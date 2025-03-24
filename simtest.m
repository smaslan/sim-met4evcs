%====================================================================================================
% This script simulates EV AC charging waveform with defined parameters (power, energy, ...) 
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
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


% show debug plots? ('': none, 'plotyy': u,i two axis plot, 'plot': u,i two separate plots, 'subplot': i,u to subplots)
cfg.dbg_plot = '';

% nominal grid voltage (rms) [V]
cfg.U_rms = 230.0;

% voltage THD [%]
cfg.U_thd_harms = 10;
cfg.U_thd = 15.0;
cfg.U_thd_mode = 'sqr';

% nominal current (rms) [A]
cfg.I_rms = 45.0;

% current THD [%]
cfg.I_thd_harms = 10;
cfg.I_thd = 10.0;
cfg.I_thd_mode = 'exp';

% PFC spurs simulator
cfg.pfc_max_h_f = 500e3;
cfg.pfc_max_h = 5;
cfg.pfc_min_f = 2e3;
cfg.pfc_max_f = 8e3;
cfg.pfc_spur_amp_rel = 0.01;

% ramp up time for current [s]
cfg.I_ramp_up_time = 0.2;
% initial current (rms) [A]
cfg.I_initial = 0.1;
% initial current time [s]
cfg.I_initial_time = 0.1;
% extra delay between current start from voltage start [s]
cfg.I_U_delay = 0.1;
% ramp down time for current [s]
cfg.I_ramp_down_time = 0.1;
% end of current to end of voltage delay [s]
cfg.I_U_end_delay = 0.1;


% steady state power factor [-]
cfg.pf = 0.9;
% or alternatively current phase shift [deg]
%cfg.I_phi = 0;

% grid frequency [Hz]
cfg.f_nom = 50.3;
% set stop frequency to generate frequency drift along the simulated interval [Hz]
cfg.f_stop = 49.9;

% digitizer sampling rate [Hz]
cfg.fs = 100000.0;

% zero padding before and after simulated waveform [s]
cfg.pad_enable = 1;
cfg.pad_init = 0.1;
cfg.pad_post = 0.1;

% simulation time (excluding padding) [s]
cfg.sim_time = 6.0;

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
[cfg.u_adc.f, cfg.u_adc.gain, cfg.u_adc.phi] = gen_adc_tfer(cfg.fs/2, 1000, u_adc.cfg, 0);

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
[cfg.i_adc.f, cfg.i_adc.gain, cfg.i_adc.phi] = gen_adc_tfer(cfg.fs/2, 1000, i_adc.cfg, 0);

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
[cfg.u_tr.f, cfg.u_tr.gain, cfg.u_tr.phi] = gen_rvd_tfer(50000.0,1000,cfg,0);

% make some current transducer transfer
cfg.Rs = 0.1;
cfg.u_Rs = cfg.Rs*10e-6;
cfg.Ls = 50e-9;
cfg.u_Ls = 5e-9;
cfg.Cp = 100e-12;
cfg.u_Cp = 10e-12;
[cfg.i_tr.f, cfg.i_tr.gain, cfg.i_tr.phi] = gen_shunt_tfer(50000.0,1000,cfg,0);

% enable ADC simulation
cfg.adc_enable = 1;

% enable transducer simulation
cfg.tr_enable = 1;

% ADC and transducer model randomization by its uncertainty?
cfg.rand_model = 0;

% generate waves
tic
[t,u,i,E_ref] = sim_evcs(cfg);
toc

% calculate energy from simulated waveforms
tic
[E_meas] = calc_energy(t,u,i,cfg);
toc

% print results
E_ref = E_ref/3600;
E_meas = E_meas/3600;
[~,str,fmt,sip,si_scale] = num_fmt(E_ref, 1e-6, 1e-6, 10, 100);
fmt = [fmt ' ' sip 'Wh'];
fprintf(['E_ref  = ' fmt '\n'],E_ref*si_scale);
fprintf(['E_meas = ' fmt '\n'],E_meas*si_scale);
fprintf(['E_dev  = ' fmt '\n'],(E_meas - E_ref)*si_scale);
rel_dev = E_ref/E_meas - 1


% figure
% x = i;
% step = ceil(0.002*adc.fs);
% window = ceil(0.003*adc.fs);
% [S,f,t] = specgram(x, 2^nextpow2(window), adc.fs, window, step);
% A = 20*log10(abs(S));
% imagesc(t,f,A);
% set(gca, 'ydir', 'normal');


 