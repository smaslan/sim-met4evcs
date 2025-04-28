%====================================================================================================
% Loader of grid impedance measurement from A1.1.3 data.
% Interpolates measured impedance to lower spots count to make processing a bit faster.
% The data are used to calculate supraharmonics voltage (or current) from current (or voltage) input.
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% Source: https://github.com/smaslan/sim-met4evcs
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
spect_file = 'Z_Grid_Urban_50th_percentile_RBW-1Hz.mat';
spect_path = fullfile(mpth,'data',spect_file);
load(spect_path);

% destination file name
model_file = strrep(spect_file,'.mat','_model.mat');
model_path = fullfile(mpth,'data',model_file);

% make complex grid Z
Z = prctile_50th_Z_Module_Ohm.*exp(j*prctile_50th_Z_Phase_rad);

% interpolate to lower spots count
f_int = linspace(0,max(f),1000);
Z_int = interp1(f, Z, f_int, 'linear', 'extrap');

% store model
model.source_file = spect_file;
model.f = f_int;
model.Z = Z_int;
save('-v6',model_path,'model');

% plot modulus
figure
plot(f,abs(Z))
hold on;
semilogx(f_int,abs(Z_int),'r')
hold off;
xlabel('f [Hz]');
ylabel('Z [\Omega]');
grid on;
box on;

% plot angle
figure
plot(f,angle(Z))
hold on;
semilogx(f_int,angle(Z_int),'r')
hold off;
xlabel('f [Hz]');
ylabel('\Phi [rad]');
grid on;
box on;



