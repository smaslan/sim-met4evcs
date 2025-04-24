%====================================================================================================
% This is demo script for function sim_evcs() that simulates EV AC charging waveform 
% with defined parameters (power, energy, ...) 
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2024 - 2025, Stanislav Maslan (smaslan@cmi.gov.cz)
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

% test mode:
%  'single' - generate waves in single slice
%  'slice' - generate waves from sub-slices of defined length
%  'both' - generate using both methods and compare difference (just for debuf purposes)
test_mode = 'single';
% slice length in samples count (for slice mode only)
slice_test_N_count = 56789;


% show debug plots? ('': none, 'plotyy': u,i two axis plot, 'plot': u,i two separate plots, 'subplot': i,u to subplots)
cfg.dbg_plot = '';

% phase count
cfg.phase_N = 2;

% optional phase shifts of voltage for each phase (default step is 360/cfg.phase_N) [deg]
cfg.U_phi = [0 120 240];

% grid frequency [Hz]
cfg.f_nom = 50.3;
% set stop frequency to generate frequency drift along the simulated interval [Hz] or leave empty or set to 0 to disable
cfg.f_stop = 49.9;

% digitizer sampling rate [Hz]
cfg.fs = 50000.0;

% nominal grid voltage (rms) [V]
%  note: scalar to set identical for all phases or vector to specify for each phase
cfg.U_rms = [230.0];

% voltage THD [%]
cfg.U_thd_harms = 10;
cfg.U_thd = 15.0;
cfg.U_thd_mode = 'sqr';

% nominal current (rms) [A]
%  note: scalar to set identical for all phases or vector to specify for each phase
cfg.I_rms = 45.0;

% current THD [%]
cfg.I_thd_harms = 10;
cfg.I_thd = 10.0;
cfg.I_thd_mode = 'exp';

% steady state power factor [-]
cfg.pf = 0.95;
% or alternatively current phase shift [deg]
%cfg.I_phi = 0;

% PFC spurs simulator
cfg.pfc_enable = 0;
cfg.pfc_max_h_f = 500e3;
cfg.pfc_max_h = 5;
cfg.pfc_min_f = 2e3;
cfg.pfc_max_f = 8e3;
cfg.pfc_spur_amp_rel = 0.01;

% initial current time interval [s]
cfg.I_initial_time = 0.1;
% initial current (rms) [A]
cfg.I_initial = 1.0;
% ramp up time for current [s]
cfg.I_ramp_up_time = 0.2;
% extra delay between current start from voltage start [s]
cfg.I_U_delay = 0.1;
% ramp down time for current [s]
cfg.I_ramp_down_time = 0.1;
% end of current to end of voltage delay [s]
cfg.I_U_end_delay = 0.1;

% simulation time (excluding padding) [s]
cfg.sim_time = 3.0;

% zero padding before and after simulated waveform [s]
cfg.pad_enable = 1;
cfg.pad_init = 0.1;
cfg.pad_post = 0.1;

% --- supra harmonics stuff:
% state
cfg.supra_enable = 1;
% lower freq limit [Hz]
cfg.supra_fmin = 1.5e3;
% upper freq limit [Hz]
cfg.supra_fmax = 150e3;
% optional generation of supraharmonics per slices (each slice will have individualy randomized frequencies/amplitudes to simulate spread spectrum/fluctuating spur)
%  define either samples count or time, or nothing to generate in one slice for entire simulation 
%cfg.supra_slice_N = 12345;
cfg.supra_slice_t = 0.1;
% optional generation of multiple spurs per time slice (alternative to generation by slices), default = 1
cfg.supra_multi_spur = 10;
% voltage emisions model path 
cfg.supra_model = fullfile(mpth, 'data', 'EV2_H2-Syn_model.mat');
% grid impedance model (for calculation of current spectrum)
cfg.supra_imp_model = fullfile(mpth, 'data', 'Z_Grid_Urban_50th_percentile_RBW-1Hz_model.mat');
% alternative grid impedance estimate if no model file provided (for calculation of current spectrum): (Rs + Ls) || Cp  
%cfg.supra_imp_Rs = 0.1;
%cfg.supra_imp_Ls = 1e-6;
%cfg.supra_imp_Cp = 1e-6;



if strcmpi(test_mode,'single')
    % optional time slice (or sample count slice) of the total simulation time (do not assign to simulate full wave at once):
    %  note: use either time values or sample count values, not both and not combination of time and samples count
    %  start of slice [s]
    %cfg.slice_t_start = 1.0;
    %  duration of slice [s] 
    %cfg.slice_t_duration = 2.0;
    %  first slice sample offset (0-based: 0-start from first sample)
    %cfg.slice_N_start = 10000;
    %  samples count for the slice 
    %cfg.slice_N_count = 10000;
end  





% enable ADC simulation
cfg.adc_enable = 1;
% enable transducer simulation
cfg.tr_enable = 1;

% ADC and transducer model randomization by its uncertainty?
cfg.rand_model = 0;

% resolution of filter being used for frequency dependent gain/phase corrections (must be x^2)
cfg.filter_size = 2^16;

% make some voltage channel ADC transfer
u_adc.cfg.dc_gain = 0.95;
u_adc.cfg.u_dc_gain = 10e-6;
u_adc.cfg.gain_fm = -500e-6;
u_adc.cfg.u_gain_fm = 500e-6;
u_adc.cfg.gain_pow = 1.5;
u_adc.cfg.gain_r_per = 20e3;
u_adc.cfg.gain_r_amp = 200e-6;
u_adc.cfg.phi_fm = -1e-3;
u_adc.cfg.u_phi_fm = 50e-6;
u_adc.cfg.u_phi_min = 2e-6;
u_adc.cfg.phi_pow = 1.0;
[cfg.u_adc.f, cfg.u_adc.gain, cfg.u_adc.phi] = gen_adc_tfer(cfg.fs/2, 100, u_adc.cfg, 0);

% make some current channel ADC transfer
i_adc.cfg.dc_gain = 1.02;
i_adc.cfg.u_dc_gain = 10e-6;
i_adc.cfg.gain_fm = -+00e-6;
i_adc.cfg.u_gain_fm = 500e-6;
i_adc.cfg.gain_pow = 1.5;
i_adc.cfg.gain_r_per = 20e3;
i_adc.cfg.gain_r_amp = 200e-6;
i_adc.cfg.phi_fm = +1e-3;
i_adc.cfg.u_phi_fm = 50e-6;
i_adc.cfg.u_phi_min = 2e-6;
i_adc.cfg.phi_pow = 1.0;
[cfg.i_adc.f, cfg.i_adc.gain, cfg.i_adc.phi] = gen_adc_tfer(cfg.fs/2, 100, i_adc.cfg, 0);

% make some voltage transducer transfer
cfg.Rhi = 100e3;
cfg.u_Rhi = cfg.Rhi*10e-6;
cfg.Chi = 0.5e-12;
cfg.u_Chi = 0.1e-12;
cfg.Rlo = 200;
cfg.u_Rlo = cfg.Rlo*10e-6;
cfg.Clo = 100e-12;
cfg.u_Clo = 10e-12;
cfg.Llo = 20e-9;
cfg.u_Llo = 5e-9;
[cfg.u_tr.f, cfg.u_tr.gain, cfg.u_tr.phi] = gen_rvd_tfer(cfg.fs/2,100,cfg,0);

% make some current transducer transfer
cfg.Rs = 0.1;
cfg.u_Rs = cfg.Rs*10e-6;
cfg.Ls = 50e-9;
cfg.u_Ls = 5e-9;
cfg.Cp = 100e-12;
cfg.u_Cp = 10e-12;
[cfg.i_tr.f, cfg.i_tr.gain, cfg.i_tr.phi] = gen_shunt_tfer(cfg.fs/2,100,cfg,0);



if ~any(strcmpi(test_mode,{'slice','single','both'}))
    error('Unknown test_mode ''%s''!',test_mode);
end



% --- generate waves (entire simulation at once):
if ~strcmpi(test_mode,'slice')
    fprintf('Generating wave in single slice...\n');
    tic
    [t,u,i,E_ref] = sim_evcs(cfg);
    toc
end


% --- generate waves using multiple slices:
if ~strcmpi(test_mode,'single')
    fprintf('Generating wave using sub-slices of length %d samples...\n',slice_test_N_count);
    tic
    % slice length (samples)
    M = slice_test_N_count;
    % total sample count
    N = (cfg.pad_enable*(cfg.pad_init + cfg.pad_post) + cfg.sim_time)*cfg.fs;    
    % no plots
    dbg_plot = cfg.dbg_plot;
    cfg.dbg_plot = '';
    
    t_slice = [];
    u_slice = [];
    i_slice = [];
    E_ref_slice = [];
    n = 0;
    while n < N
        %  first slice sample offset (0-based: 0=start from first sample)
        cfg.slice_N_start = n;
        %  samples count for the slice 
        cfg.slice_N_count = M;
        
        % generate slice
        [tx,ux,ix,E_ref_slice(end+1,:)] = sim_evcs(cfg);
                
        % merge slices
        t_slice = [t_slice; tx];
        u_slice = [u_slice; ux];
        i_slice = [i_slice; ix];
        
        % next slice
        n = n + M;
    end
    E_ref_slice = sum(E_ref_slice,1);
    cfg.dbg_plot = dbg_plot;
    toc
    
    if strcmpi(test_mode,'slice')
        % use slice output as test result
        t = t_slice;
        u = u_slice;
        i = i_slice;
        E_ref = E_ref_slice;        
    
    elseif strcmpi(test_mode,'both')
        % compare single and slice modes
        
        dev_t_slice = sum(abs(t - t_slice))
        dev_u_slice = sum(abs(u - u_slice))
        dev_i_slice = sum(abs(i - i_slice))
        dev_E_ref = (E_ref - E_ref_slice)/3600
                
        % plot wave differences
        if strcmp(cfg.dbg_plot,'plotyy')
            leg = {};
            for phid = 1:cfg.phase_N
                leg{end+1} = sprintf('U(L%d)',phid);            
            end
            for phid = 1:cfg.phase_N
                leg{end+1} = sprintf('I(L%d)',phid);            
            end
            
            figure;
            [hax] = plotyy(t,u - u_slice, t,i - i_slice);
            xlabel(hax(1), 'time [s]');
            ylabel(hax(1), 'u - u\_slice [V]');
            ylabel(hax(2), 'i - i\_slice [A]');
            grid on;
            box on;
            legend(leg);
        
        elseif strcmp(cfg.dbg_plot,'plot')
            leg = {};
            for phid = 1:cfg.phase_N
                leg{phid} = sprintf('L%d',phid);            
            end
            
            figure;
            plot(t, u - u_slice);
            xlabel('time [s]');
            ylabel('u - u\_slice [V]');
            grid on;
            box on;
            legend(leg);
            
            figure;
            plot(t, i - i_slice);
            xlabel('time [s]');
            ylabel('i - i\_slice [A]');
            grid on;
            box on;
            legend(leg);
            
        elseif strcmp(cfg.dbg_plot,'subplot')
            leg = {};
            for phid = 1:cfg.phase_N
                leg{phid} = sprintf('L%d',phid);            
            end
        
            figure;
            subplot(2,1,1);
            plot(t, u - u_slice);
            xlabel('time [s]');
            ylabel('u - u\_slice [V]');
            grid on;
            box on;
            legend(leg);
            
            subplot(2,1,2);        
            plot(t, i - i_slice);
            xlabel('time [s]');
            ylabel('i - i\_slice [A]');
            grid on;
            box on;
            legend(leg);        
        
        end
    end
    
    clear t_slice u_slice i_slice; 
    
end



  

% --- calculate energy from simulated waveforms:
fprintf('Calculating energy of generated waves...\n');
tic
[E_meas] = calc_energy(t,u,i,cfg);
toc

% print results
E_ref = E_ref/3600;
E_meas = E_meas/3600;
E_ref_sum = sum(E_ref);
E_meas_sum = sum(E_meas);
[~,str,fmt,sip,si_scale] = num_fmt(max(E_ref), 1e-6, 1e-6, 10, 100);
fmt = [fmt ' ' sip 'Wh '];
fprintf(['E_ref  = [%s]\n'],sprintf(fmt, E_ref*si_scale));
fprintf(['E_meas = [%s]\n'],sprintf(fmt, E_meas*si_scale));
fprintf(['E_dev  = [%s]\n'],sprintf(fmt, (E_meas - E_ref)*si_scale));

if cfg.phase_N > 1
    [~,str,fmt,sip,si_scale] = num_fmt(max(E_ref_sum), 1e-6, 1e-6, 10, 100);
    fmt = [fmt ' ' sip 'Wh '];
    fprintf(['E_ref_sum  = %s\n'],sprintf(fmt, E_ref_sum*si_scale));
    fprintf(['E_meas_sum = %s\n'],sprintf(fmt, E_meas_sum*si_scale));
    fprintf(['E_dev_sum  = %s\n'],sprintf(fmt, (E_meas_sum - E_ref_sum)*si_scale));
end
rel_dev = E_ref_sum/E_meas_sum - 1


% figure
% x = i;
% step = ceil(0.002*adc.fs);
% window = ceil(0.003*adc.fs);
% [S,f,t] = specgram(x, 2^nextpow2(window), adc.fs, window, step);
% A = 20*log10(abs(S));
% imagesc(t,f,A);
% set(gca, 'ydir', 'normal');


 