function [t,u,i,E_sim,u_raw,i_raw] = sim_evcs(cfg)
% This script simulates EV AC charging waveforms with defined parameters (U, I, PF, ...)
% Returns simulated waveforms and calculated active energy. 
% It can also apply frequency dependent gain-phase model of ADC and transducers and return modified waveforms.
% Frequency dependent corrections are performed using FFT filtering. 
%
% Usage:
%   [t,u,i,E,u_raw,i_raw] = sim_evcs(cfg)
%
% Parameters:
%   cfg.dbg_plot - optional, show debug plot? {'': no, 'plotyy' two axis mode, 'plot': two plots, 'subplot': subplots}
%   cfg.U_rms - effective voltage [V]
%   cfg.I_rms - effective current [A]
%   cfg.U_phi - phase angle for both U and I [deg] (for multiphase systems)
%   cfg.I_phi or cfg.pf - current to voltage phase shift [deg] of effective power factor
%   cfg.f_nom - nominal grid frequency [Hz]
%   cfg.f_stop - optional, target frequency [Hz], will generate linear frequency sweep from f_nom to f_stop
%   cfg.U_thd_mode - optional, shape of fundamental harmonics {'exp': exponential decay, 'sqr': square wave}
%                  - following parameters are required when enabled:
%     cfg.U_thd_harms - number of harmonics to generate, e.g. 7: up to 7th harmonic 
%     cfg.U_thd - required fundamental referenced THD value [%]
%   cfg.I_thd_mode - optional, shape of fundamental harmonics {'exp': exponential decay, 'sqr': square wave}
%                  - following parameters are required when enabled:
%     cfg.I_thd_harms - number of harmonics to generate, e.g. 7: up to 7th harmonic 
%     cfg.I_thd - required fundamental referenced THD value [%]
%   cfg.pfc_enable - optional, enable generation of PFC-like harmonics
%                  - generates triangular wave with variable frequency emulating CrM mode, requires following:
%     cfg.pfc_max_h_f - upper frequency limit for harmonics [Hz]
%     cfg.pfc_max_h - maximum harmonic id of pfc switching frequency fundamental, e.g. 7: max 7th harmonic
%     cfg.pfc_min_f - lowest switching frequency [Hz]
%     cfg.pfc_max_f - highest switching frequency [Hz] 
%     cfg.pfc_spur_amp_rel - relative amplitude of spur to fundamental grid frequency component [-]
%   cfg.sim_time - simulation time excluding padding [s]
%   cfg.fs - sampling rate [Hz]
%   cfg.pad_enable - optional, waveform padding enable 
%   cfg.pad_init - optional, initial waveform padding [s] 
%   cfg.pad_post - optional, ending waveform padding [s]
%   cfg.I_ramp_up_time - optional, current wave rampup duration [s]
%   cfg.I_initial - optional, initial current (before ramp) [A]
%   cfg.I_initial_time - optional, initial current time (before ramp) [s]
%   cfg.U_I_delay - optional, extra delay between voltage start and initial current [s]
%   cfg.I_ramp_down_time - optional, ramp-down time [s]
%   cfg.I_U_end_delay - optional, extra delay between current end and voltage end [s]
%  System model related stuff:
%   cfg.filter_size - optional, FFT filter mask resolution used for the frequency dependent gain, phase corrections
%                     must be 2^x size
%   cfg.tfer_interp - optional, default interpolation mode for ADC and transducer transfers (default: 'pchip')
%   cfg.adc_enable - enable ADC transfer model (will apply cfg.u_adc and cfg.i_adc transfers to generated waveforms)
%   cfg.u_adc - gain-phase transfer of voltage channel ADC:
%     cfg.u_adc.f - frequency vector [Hz] ranging from 0 to cfg.fs/2
%     cfg.u_adc.gain.v - gain for each f [V/V]
%     cfg.u_adc.gain.u - gain absolute uncertainty for each f [V/V] 
%     cfg.u_adc.phi.v - phase for each f [rad]
%     cfg.u_adc.phi.u - phase absolute uncertainty for each f [rad]
%   cfg.i_adc - gain-phase transfer of current channel ADC, items same as cfg.u_adc   
%   cfg.tr_enable - enable transducer transfer model (will apply cfg.u_tr and cfg.i_tr transfers to generated waveforms)
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
%   cfg.rand_model - randomize ADC and transducer models by their uncertainty (useful for monte-carlo simulations)
%
% Returns:
%   cfg.t - time vector of generated waveforms [s]
%   cfg.u - voltage waveform [V]
%   cfg.i - current waveform [A]
%   cfg.E_sim - active energy of simulated u/i waveforms [Ws]
%   cfg.u_raw - raw voltage waveform before applying system model [V] 
%   cfg.i_raw - raw current waveform before applying system model [A]
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.
%

    % show debug plots? ('': none, 'plotyy': u,i two axis plot, 'plot': u,i two separate plots, 'subplot': i,u to subplots)
    cfg = def(cfg, 'dbg_plot', '');
    
    cfg = def(cfg, 'U_thd_mode', '');
    if ~isempty(cfg.U_thd_mode)
        check_str(cfg, {'U_thd_harms', 'U_thd'});        
    end
    
    cfg = def(cfg, 'I_thd_mode', '');
    if ~isempty(cfg.I_thd_mode)
        check_str(cfg, {'I_thd_harms', 'I_thd'});        
    end
    
    cfg = def(cfg, 'pfc_enable', 0);
    if cfg.pfc_enable
        check_str(cfg, {'pfc_max_h_f', 'pfc_max_h', 'pfc_min_f', 'pfc_max_f', 'pfc_spur_amp_rel'});
    end
    
    % default phase angle
    cfg = def(cfg, 'U_phi', 0);
    U_phi = cfg.U_phi*180.0/pi;
    
    % default current ramps
    cfg = def(cfg, 'I_ramp_up_time', 0.0);
    cfg = def(cfg, 'I_initial', 0.0);
    cfg = def(cfg, 'I_initial_time', 0.0);
    cfg = def(cfg, 'I_U_delay', 0.0);
    cfg = def(cfg, 'I_ramp_down_time', 0.0);
    cfg = def(cfg, 'I_U_end_delay', 0.0);
    
    % select current phase angle
    if ~xor(isfield(cfg, 'pf'), isfield(cfg, 'I_phi'))
        error('Use either power factor ''pf'' or current phase shift ''I_phi''!');
    elseif isfield(cfg, 'pf')
        I_phi = mod(acos(-cfg.pf) + pi/2, pi) - pi/2;
    elseif isfield(cfg, 'I_phi')
        I_phi = cfg.I_phi/180.0*pi;
    end
    I_phi = I_phi + U_phi;
    
    % default grid frequency [Hz]
    cfg = def(cfg, 'f_nom', 50.0);    
    
    % default padding
    cfg = def(cfg, 'pad_enable', 0);
    cfg = def(cfg, 'pad_init', 0.0);
    cfg = def(cfg, 'pad_post', 0.0);
       
    % default resolution of filter being used for frequency dependent gain/phase corrections (must be x^2)
    cfg = def(cfg, 'filter_size', 8192);
    
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
    
    % total simulated time
    t_sim = cfg.sim_time;
    t_start = 0.0;
    t_stop = cfg.sim_time;
    if cfg.pad_enable
        t_start = cfg.pad_init;
        t_stop = t_start + t_sim;
        t_sim = t_sim + cfg.pad_init + cfg.pad_post;
    end
    
    % time vector
    t = [];
    t(:,1) = [0:1/cfg.fs:t_sim-1/cfg.fs];
    
    % generate fund. frequency vector along simulated interval
    w0 = 2*pi*cfg.f_nom;
    if isfield(cfg,'f_stop') && cfg.f_stop
        w0 = 2*pi*linspace(cfg.f_nom, cfg.f_stop, numel(t)).';
    end
    
    % fundamental grid voltage component
    A = cfg.U_rms*2^0.5;
    
    % generate grid voltage harmonics 
    [f_uh, A_uh, ph_uh] = gen_thd_harms(cfg.U_thd, cfg.U_thd_mode, cfg.U_thd_harms);
    A_uh = A*A_uh;
    
    % fix voltage components to match desired rms
    ku_rms = cfg.U_rms/rms(2^-0.5*[A A_uh]);
    A = A*ku_rms;
    A_uh = A_uh*ku_rms;
    
    % list of all voltage harmonics
    f_list = [1 f_uh];
    A_list = [A A_uh];
    ph_list = [0 ph_uh] + U_phi;
    
    
    
    % generate voltage harmonics (one by one to save memory for large datasets)
    u = zeros(size(t));
    for k = 1:numel(f_list)
        wt = w0.*f_list(k).*(t - t_start);
        u = u + A_list(k)*sin(wt + ph_list(k));
    end
    clear wt;
    
    % apply voltage padding
    u(find(t < t_start | t > t_stop)) = 0.0;
    
    
    
    % make current envelope
    t_i_start = t_start + cfg.I_U_delay;
    t_i_end_initial = t_i_start + cfg.I_initial_time;
    t_i_end_rampup = t_i_end_initial + cfg.I_ramp_up_time;
    t_i_end = t_stop - cfg.I_U_end_delay;
    t_i_end_nom = t_i_end - cfg.I_ramp_down_time;
    t_i_steps = [0, t_i_start, t_i_start, t_i_end_initial, t_i_end_rampup, t_i_end_nom, t_i_end, t_sim];
    i_env_steps = [0, 0, cfg.I_initial, cfg.I_initial, cfg.I_rms, cfg.I_rms, 0, 0]/cfg.I_rms;
    t_i_steps = t_i_steps + eps*[1:numel(t_i_steps)];
    i_env = interp1(t_i_steps, i_env_steps, t,'linear','extrap');
    
    
    
    % main current component
    A = cfg.I_rms*2^0.5;
     
     % generate current THD components
    [f_uh, A_uh, ph_uh] = gen_thd_harms(cfg.I_thd, cfg.I_thd_mode, cfg.I_thd_harms);
    A_uh = A*A_uh;
    
    % fix current components to match desired rms
    ki_rms = cfg.I_rms/rms(2^-0.5*[A A_uh]);
    A = A*ki_rms;
    A_uh = A_uh*ki_rms;
    
    % make list of current harmonics
    f_list = [1 f_uh];
    A_list = [A A_uh];
    ph_list = [0 ph_uh] + I_phi;
    
    % generate current harmonics (one by one to save memory for large datasets)
    i = zeros(size(t));
    for k = 1:numel(f_list)
        wt = w0.*f_list(k).*(t - t_start);
        i = i + A_list(k)*sin(wt + ph_list(k)).*i_env;
    end
    clear wt;
    
    % generate PFC spurs
    if cfg.pfc_enable
        ip = gen_pfc_emi(t - t_start, w0.*(t - t_start), cfg.pfc_min_f,cfg.pfc_max_f, cfg.pfc_max_h,cfg.pfc_max_h_f, A.*i_env*cfg.pfc_spur_amp_rel);
        i = i + ip;
    end    
    clear i_env;
        
    % calculate actual simulated energy dose [Ws]
    p = u.*i;
    E_sim = mean(p)*t_sim;
    
    % return unmodified waveforms
    u_raw = u;
    i_raw = i;
      
    
    
    
    % --- Apply ADC and transducer models:
    
    % randomize system model by uncertainty?
    cfg = def(cfg, 'rand_model', 0);
    cfg.rand_model = ~~cfg.rand_model;      
    
    % default ADC/transducer transfer interpolation mode
    cfg = def(cfg, 'tfer_interp', 'pchip');
        
    if cfg.adc_enable || cfg.tr_enable
    
        % define channels to process 
        channels{1}.y = u; % input
        channels{1}.y_var = 'u'; % output variable
        channels{1}.adc_tfer = cfg.u_adc; % adc model
        channels{1}.tr_tfer = cfg.u_tr; % transducer model
        channels{2}.y = i; % input
        channels{2}.y_var = 'i'; % output variable
        channels{2}.adc_tfer = cfg.i_adc; % adc model
        channels{2}.tr_tfer = cfg.i_tr; % transducer model
        
        % for each channel:
        for k = 1:numel(channels)
            chn = channels{k};
        
            % get combined frequency vector of ADC and transducer correction
            tf_freq = [];
            if cfg.adc_enable
                tf_freq = chn.adc_tfer.f(:);
            end
            if cfg.tr_enable
                tf_freq = [tf_freq;chn.tr_tfer.f(:)];
            end
            tf_freq = unique(sort(tf_freq));
            if min(tf_freq) > 0 || max(tf_freq) < cfg.fs/2
                error('Frequency vector of ADC and transducer transfers must cover frequency range from 0 to nyquist (cfg.fs/2)!');
            end
            
            % combine ADC and transducer transfers:
            %   note: optional randomization of tfers assuming correlated uncertainty along frequency axis            
            
            % ADC transfer
            tf_gain = ones(size(tf_freq));
            tf_phi = zeros(size(tf_freq));
            if cfg.adc_enable
                freq = chn.adc_tfer.f(:);
                gain = chn.adc_tfer.gain.v(:) + cfg.rand_model*randn()*chn.adc_tfer.gain.u(:);    
                phi = chn.adc_tfer.phi.v(:) + cfg.rand_model*randn()*chn.adc_tfer.phi.u(:);            
                tf_gain = interp1(freq, gain, tf_freq, cfg.tfer_interp, 'extrap');
                tf_phi = interp1(freq, phi, tf_freq, cfg.tfer_interp, 'extrap');
            end
            
            % transducer transfer
            if cfg.tr_enable
                freq = chn.tr_tfer.f(:);
                gain = chn.tr_tfer.gain.v(:) + cfg.rand_model*randn()*chn.tr_tfer.gain.u(:);    
                phi = chn.tr_tfer.phi.v(:) + cfg.rand_model*randn()*chn.tr_tfer.phi.u(:);            
                gain = interp1(freq, gain, tf_freq, cfg.tfer_interp, 'extrap');
                phi = interp1(freq, phi, tf_freq, cfg.tfer_interp, 'extrap');
                tf_gain = tf_gain.*gain;
                tf_phi = tf_phi + phi;
            end
            
            % add extra pading needed for filter
            filter_pad = cfg.filter_size;
            pad = zeros(filter_pad,1);
            y_pad = [pad;chn.y;pad];
            
            % apply filter
            [y_filt, id_start, id_stop] = td_fft_filter(y_pad, cfg.fs, cfg.filter_size, tf_freq, tf_gain, tf_phi);
                        
            % get rid of padding residue
            eval([chn.y_var ' = y_filt((filter_pad - id_start + 1):(filter_pad - id_start + numel(t)));']);
                
            clear y_pad y_filt;
        end
    
    end 
    
    
    % show some debug plots?
    if strcmp(cfg.dbg_plot,'plotyy')
        figure;
        [hax] = plotyy(t,u_raw, t,i_raw);
        xlabel(hax(1), 'time [s]');
        ylabel(hax(1), 'u [V]');
        ylabel(hax(2), 'i [A]');
        grid on;
        box on;
    
    elseif strcmp(cfg.dbg_plot,'plot')
        figure;
        plot(t, u_raw);
        xlabel('time [s]');
        ylabel('u [V]');
        grid on;
        box on;
        
        figure;
        plot(t, i_raw);
        xlabel('time [s]');
        ylabel('i [A]');
        grid on;
        box on;
        
    elseif strcmp(cfg.dbg_plot,'subplot')
        figure;
        subplot(2,1,1);
        plot(t, u_raw);
        xlabel('time [s]');
        ylabel('u [V]');
        grid on;
        box on;
        
        subplot(2,1,2);        
        plot(t, i_raw);
        xlabel('time [s]');
        ylabel('i [A]');
        grid on;
        box on;        
    
    end 
    
    
    % figure
    % x = i;
    % step = ceil(0.002*cfg.fs);
    % window = ceil(0.003*cfg.fs);
    % [S,f,t] = specgram(x, 2^nextpow2(window), cfg.fs, window, step);
    % A = 20*log10(abs(S));
    % imagesc(t,f,A);
    % set(gca, 'ydir', 'normal');


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