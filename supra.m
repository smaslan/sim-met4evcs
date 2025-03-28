%====================================================================================================
% Supraharmonics model generator from data from A1.1.3 wave records.
% Processes raw waveform, identifies noise level and significant spurs.
% Generates model describing the noise and spurs.
% Reconstructs the signal from model. 
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2025, Stanislav Maslan (smaslan@cmi.gov.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.
%====================================================================================================
clc;
clear all;
close all;

% current path
mpth = fileparts(mfilename('fullpath'));
cd(mpth);
addpath(mpth);

if isOctave
    % set QT as default Octave graphs toolkit
    % ###note: for some reason, in Windows it is available only if Octave started with octave-gui.exe and not by octave-cli.exe
    graphics_toolkit('qt');
end

spec_file_folder = 'data';

% source file
spect_file = 'EV2_H2-Syn.mat';
spect_path = fullfile(mpth,'data',spect_file);
load(spect_path);

% destination file name
model_file = strrep(spect_file,'.mat','_model.mat');
model_path = fullfile(mpth,'data',model_file);

% limit max processed frequency [Hz]
f_max = 500e3;

% split FFT analysis to time slices (windows)
%  note: calculating spectra from small slices will give more information of peak amplitudes compared to spectrum of entire record, where
%        the spead spectrum components will have lower average amplitudes
N_slice = 200e3;

% split wave to multiple time slices
N = numel(A);
S = floor(N/N_slice);
A_slices = A(1:S*N_slice);
A_slices = reshape(A_slices, [N_slice,S]);

% calculate amplitude spectra of slices
fs = 1/Tinterval;
N = N_slice;
amp = fft(A_slices);
fax(:,1) = [0:floor(N/2)-1]/N*fs;
amp = amp(1:floor(N/2),:)*2/N;
amp = amp(fax < f_max,:);
fax = fax(fax < f_max);
amp = abs(amp);
% find amp maxima of all slices
amp = max(amp,[],2);
N = numel(amp);

% find maximum and minimum levels from sub-slices
spur_filt_bw = 500;
MF = round(spur_filt_bw/diff(fax(1:2)));
S = floor(N/MF);
amp_m = amp(1:S*MF);
amp_m = reshape(amp_m,[MF,S]);
amp_max = max(amp_m,[],1)';
amp_min = min(amp_m,[],1)';
amp_max_f = fax(round(MF/2):MF:end);
amp_max_f = amp_max_f(1:S);
amp_max_f = [min(fax);amp_max_f;max(fax)];
amp_max = [amp_max(1);amp_max;amp_max(end)];
amp_min = [amp_min(1);amp_min;amp_min(end)];

% find median levels
median_slice_bw = 5000;
MF = round(median_slice_bw/diff(fax(1:2)));;
S = floor(N/MF);
amp_m = amp(1:S*MF);
amp_m = reshape(amp_m,[MF,S]);
amp_med = median(amp_m,1)';
amp_med_f = fax(MF/2:MF:end);
amp_med_f = amp_med_f(1:S);
amp_med_f = [min(fax);amp_med_f;max(fax)];
amp_med = [amp_med(1);amp_med;amp_med(end)];

% filter min/max levels from sub-slices to get estimate of noise (excluding peaks)
filt_bw = 15000;
MF = round(filt_bw/diff(amp_max_f(2:3)));
P = numel(amp_max);
S = floor(P/MF);
amp_m = amp_max(1:S*MF);
amp_m = reshape(amp_m,[MF,S]);
amp_med_max = median(amp_m,1)';
amp_med_max_f = amp_max_f(round(MF/2):MF:end);
amp_med_max_f = amp_med_max_f(1:S);
amp_med_max_f = [min(fax);amp_med_max_f;max(fax)];
amp_med_max = [amp_med_max(1);amp_med_max;amp_med_max(end)];
amp_m = amp_min(1:S*MF);
amp_m = reshape(amp_m,[MF,S]);
amp_med_min = median(amp_m,1)';
amp_med_min = [amp_med_min(1);amp_med_min;amp_med_min(end)];

% detection threshold for peaks over noise
spur_amp_thr = interp1(amp_med_max_f,amp_med_max, amp_max_f, 'linear', 'extrap')*1.15;
sid = [1:numel(amp_max)];
sid = sid(amp_max > spur_amp_thr);
spur_N = numel(sid)
% extract spurs
spur_fx = amp_max_f(sid);
spur_ax = amp_max(sid);
spur_ufx = 0.5*spur_filt_bw;

% store model
model.source_file = spect_file;
model.noise_f = amp_med_max_f;
model.noise_a_max = amp_med_max;
model.noise_a_min = amp_med_max;
model.noise_a_med = interp1(amp_med_f,amp_med, amp_med_max_f, 'linear', 'extrap');
model.spur_f = spur_fx;
model.spur_uf = spur_ufx;
model.spur_a = spur_ax;
save('-v6',model_path,'model');

 

% plot input analysis
figure(1);
semilogy(fax, amp)
hold on;
semilogy(amp_max_f, amp_max, 'r')
semilogy(amp_max_f, amp_min, 'r')
semilogy(amp_med_f, amp_med, 'm')
semilogy(amp_med_max_f, amp_med_max, 'k')
semilogy(amp_med_max_f, amp_med_min, 'k')
semilogy(amp_max_f, spur_amp_thr, 'g')
semilogy(amp_max_f(sid), amp_max(sid), 'ro')
hold off;
title('Source waveform');
xlabel('f [Hz]');
ylabel('A [-]');
grid on;
box on;


% --- reconstruct signal from just obtained model:

% model sample count
N_sim = 1000000;
N = ceil(N_sim/2)*2; % always even sample count, at the end will be reduced to desired count
% model sampling rate [Hz]
fs_sim = 2e6;
% model freq limits [Hz]
f_sim_max = 500e3;
f_sim_min = 2e3;

% model freq axis
f_sim(:,1) = [0:N/2]/N*fs_sim;
f_step = diff(f_sim(1:2));
f_mask = (f_sim >= f_sim_min & f_sim <= f_sim_max);
f_mask(1) = 0; % no DC spot
f_mask(end) = 0; % no nyquist spot

% generate noise level       
noise_f = f_sim;
noise_amp = interp1(amp_med_f,amp_med, noise_f, 'linear', 'extrap');
noise_max = interp1(amp_med_max_f,amp_med_max, noise_f, 'linear', 'extrap');
noise_min = interp1(amp_med_max_f,amp_med_min, noise_f, 'linear', 'extrap');
noise_std = (noise_max - noise_amp);
amp_sim = f_mask.*max(noise_amp + noise_std.*2.*(rand(size(noise_f)) - 0.5),noise_min);

% optional multiple spurs to simulate spread spectrum
multi_spur = 1;
% now generate spurs with randomized frequencies corresponding to analysis freq resolution
for k = 1:multi_spur
    spur_bin = round(rand2(size(spur_fx),spur_fx,spur_ufx)/f_step);
    amp_sim(spur_bin) = spur_ax;
end

% figure(2);
% semilogy(fax, amp)
% hold on;
% semilogy(noise_f, amp_sim, 'r')
% hold off;

% mirror simulated spectrum to negative frequencies
amp_sim = [amp_sim;conj(amp_sim(end-1:-1:2))];
sim_wave = real(ifft(amp_sim))*numel(amp_sim)/2;
sim_wave = sim_wave(1:N_sim);

% rerun spectrum analysis for validation
amp_mod(:,1) = fft(sim_wave);
N_mod = numel(amp_mod);
amp_mod = amp_mod(1:floor(N_mod/2))*2/N_mod;
fax_mod(:,1) = [0:floor(N_mod/2)-1]/N_mod*fs_sim;
amp_mod = amp_mod(fax_mod < f_max & fax_mod > f_sim_min);
fax_mod = fax_mod(fax_mod < f_max & fax_mod > f_sim_min);
amp_mod = abs(amp_mod);


figure(3);
semilogy(fax_mod, amp_mod)
title('Simulated waveform');
xlabel('f [Hz]');
ylabel('A [-]');
grid on;
box on;

