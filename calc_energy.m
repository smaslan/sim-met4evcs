function [E] = calc_energy(t,u,i,cfg)
% Simple algorithm for measuring energy of waveform.
% Also applies ADC and transducer corrections.
% Note: it is just simple testing algorithm for sim_evcs() testing.
%
% Usage:
%   [E] = calc_energy(t,u,i,cfg)
%
% Parameters:
%   t - sampling time vector [s]
%   u - digitized voltage waveform [V]
%   i - digitizer current waveform [V]
%   cfg.filter_size - optional, FFT filter mask resolution used for the frequency dependent gain, phase corrections
%   cfg.adc_enable - enable ADC correction:
%   cfg.u_adc - gain-phase transfer of voltage channel ADC:
%     cfg.u_adc.f - frequency vector [Hz] ranging from 0 to cfg.fs/2
%     cfg.u_adc.gain.v - gain for each f [V/V]
%     cfg.u_adc.gain.u - gain absolute uncertainty for each f [V/V] 
%     cfg.u_adc.phi.v - phase for each f [rad]
%     cfg.u_adc.phi.u - phase absolute uncertainty for each f [rad]
%   cfg.i_adc - gain-phase transfer of current channel ADC, items same as cfg.u_adc   
%   cfg.tr_enable - enable transducer transfer model:
%   cfg.u_tr - gain-phase transfer of voltage channel transducer:
%     cfg.u_tr.f - frequency vector [Hz] ranging from 0 to cfg.fs/2
%     cfg.u_tr.gain.v - output/input gain for each f [V/V]
%     cfg.u_tr.gain.u - output/input gain absolute uncertainty for each f [V/V] 
%     cfg.u_tr.phi.v - output to input phase for each f [rad]
%     cfg.u_tr.phi.u - output to input phase absolute uncertainty for each f [rad]
%   cfg.i_tr - gain-phase transfer of current channel transducer:
%     cfg.i_tr.f - frequency vector [Hz] ranging from 0 to cfg.fs/2
%     cfg.i_tr.gain.v - shunt impedance modulus |Z| for each f [Ohm]
%     cfg.i_tr.gain.u - shunt impedance modulus absolute uncertainty for each f [Ohm] 
%     cfg.i_tr.phi.v - shunt impedance phase angle for each f [rad]
%     cfg.i_tr.phi.u - shunt impedance phase angle absolute uncertainty for each f [rad]
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.
%

    % sampling rate [Hz] assuming equidistant sampling
    fs = 1/(t(2) - t(1));
    
    % total sampling time [s]
    t_meas = t(end) - t(1) + 1/fs;
    
    % default ADC transfers
    cfg = def(cfg, 'adc_enable', 0);
    if cfg.adc_enable
        check_str(cfg, {'u_adc', 'i_adc'});
        check_str(cfg.u_adc, {'f', 'gain', 'phi'});
        check_str(cfg.i_adc, {'f', 'gain', 'phi'});
    else
        cfg.u_adc = [];
        cfg.i_adc = [];
    end    
    
    % default transducer transfers
    cfg = def(cfg, 'tr_enable', 0);
    if cfg.tr_enable
        check_str(cfg, {'u_tr', 'i_tr'});
        check_str(cfg.u_tr, {'f', 'gain', 'phi'});
        check_str(cfg.i_tr, {'f', 'gain', 'phi'});
    else
        cfg.u_tr = [];
        cfg.i_tr = [];
    end
    
    % default resolution of filter being used for frequency dependent gain/phase corrections (must be x^2)
    cfg = def(cfg, 'filter_size', 8192);
    if abs(rem(log2(cfg.filter_size),1)) > eps
        error('filter_size must be of size 2^x!');
    end
    
    % add extra pading needed for filter
    filter_pad = cfg.filter_size;
    pad = zeros(filter_pad,1);
    
    if cfg.adc_enable
        
        % voltage ADC transfer 
        tf_freq = cfg.u_adc.f;
        tf_gain = 1./cfg.u_adc.gain.v;
        tf_phi = -cfg.u_adc.phi.v;
        y_pad = [pad;u;pad];
        [u, id_start, id_stop] = td_fft_filter(y_pad, fs, cfg.filter_size, tf_freq, tf_gain, tf_phi);
        u = u((filter_pad - id_start + 1):(filter_pad - id_start + numel(t)));            
        
        % current ADC transfer 
        tf_freq = cfg.i_adc.f;
        tf_gain = 1./cfg.i_adc.gain.v;
        tf_phi = -cfg.i_adc.phi.v;
        y_pad = [pad;i;pad];
        [i, id_start, id_stop] = td_fft_filter(y_pad, fs, cfg.filter_size, tf_freq, tf_gain, tf_phi);
        i = i((filter_pad - id_start + 1):(filter_pad - id_start + numel(t)));
        
    end      
    
    if cfg.tr_enable
        
        % voltage transducer transfer 
        tf_freq = cfg.u_tr.f;
        tf_gain = 1./cfg.u_tr.gain.v;
        tf_phi = -cfg.u_tr.phi.v;
        y_pad = [pad;u;pad];
        [u, id_start, id_stop] = td_fft_filter(y_pad, fs, cfg.filter_size, tf_freq, tf_gain, tf_phi);
        u = u((filter_pad - id_start + 1):(filter_pad - id_start + numel(t)));            
        
        % current transducer transfer 
        tf_freq = cfg.i_tr.f;
        tf_gain = 1./cfg.i_tr.gain.v;
        tf_phi = -cfg.i_tr.phi.v;
        y_pad = [pad;i;pad];
        [i, id_start, id_stop] = td_fft_filter(y_pad, fs, cfg.filter_size, tf_freq, tf_gain, tf_phi);
        i = i((filter_pad - id_start + 1):(filter_pad - id_start + numel(t)));
        
    end
    
    % calculate actual simulated energy dose [Ws]
    p = u.*i;
    E = mean(p)*t_meas;

end


function str = def(str, item, value)
% create default structure item value if not exist
    if ~isfield(str, item)
        str = setfield(str, item, value);
    end
end

function [] = check_str(str, list)
% check missing structure items, generate error if missing some
    miss = ~isfield(str, list);
    if any(miss)        
        error(sprintf('Missing parameters: %s',catcellcsv(list(miss),', ')));
    end    
end 