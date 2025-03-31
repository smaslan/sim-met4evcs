function [f,qu_gain,qu_phi] = gen_adc_tfer(f_max,f_count, cfg, varargin)
% Simple function to simulate parametric frequency transfer of digitizer gain and phase.
%
% Usage:
%   [f,qu_gain,qu_phi] = gen_adc_tfer(f_max,f_count, cfg) 
%   [f,qu_gain,qu_phi] = gen_adc_tfer(f_max,f_count, cfg, debug_plot)
%
% Parameters:
%   f_max - max frequency of freq transfer (usually nyquist fs/2)
%   f_count - number of frequency spots to generate
%   cfg - transfer parametrs:
%   cfg.dc_gain - absolute gain at DC [V/V]
%   cfg.u_dc_gain - DC gain absolute uncertainty [V/V]
%   cfg.gain_fm - gain deviation from 1.0 at f_max [V/V]
%   cfg.gain_pow - frequency factor of gain error (1: linear, 2: with square of f, etc.)
%   cfg.u_gain_fm - absolute gain uncertainty at f_max [V/V]
%   cfg.gain_r_per - optional gain ripple period [Hz] (simulates FIR filter ripple)
%   cfg.gain_r_amp - optional gain ripple amplitude [V/V]
%   cfg.phi_fm - phase error at f_max [rad]
%   cfg.phi_pow - frequency factor of phase error (1: linear, 2: with square of f, etc.)
%   cfg.u_phi_fm - absolute phase uncertainty at f_max [rad]
%   cfg.u_phi_min - absolute phase uncertainty at f near 0 [rad]
%   debug_plot - optional plot of generated tfers
%
% Returns:
%   f - generated frequency vector from 0 to f_max [Hz]
%   qu_gain.v - gain transfer [V/V]   
%   qu_gain.u - absolute gain transfer uncertainty [V/V]
%   qu_phi.v - phase transfer [rad]   
%   qu_phi.u - absolute phase transfer uncertainty [rad]
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.

    % debug plot?
    do_plot = nargin() > 3 && isnumeric(varargin{1}) && varargin{1};

    % generate freq. vector
    f_low = 1000.0;
    f = [0, logspace(log10(1),log10(min(f_low,f_max)), round(f_count/2)), linspace(min(f_low+1, f_max),f_max + eps,round(f_count/2))].';
    f = unique(sort(f));
    
    % generate smooth gain tfer   
    gain(:,1) = cfg.dc_gain*(1 + cfg.gain_fm.*(f/max(f)).^cfg.gain_pow);
    
    % generate uncertainty of gain tfer
    u_gain(:,1) = (cfg.u_dc_gain^2 + (cfg.u_gain_fm.*(f/max(f)).^cfg.gain_pow).^2).^0.5;
    
    if isfield(cfg,'gain_r_per') && cfg.gain_r_per && isfield(cfg,'gain_r_amp') && cfg.gain_r_amp
        % generate gain tfer ripple (5922 FIR filter-like shape)
        gain_osc = cfg.gain_r_amp.*sin(f/cfg.gain_r_per*2*pi);
        % add to smooth gain
        gain = gain + gain_osc;
    end
    
    % return gain tfer quantity
    qu_gain.v = gain;
    qu_gain.u = u_gain;
        
    if do_plot
        figure;
        semilogx(f,gain)
        hold on;
        semilogx(f,gain + u_gain,'r')
        semilogx(f,gain - u_gain,'r')
        xlabel('f [Hz]')
        ylabel('gain [V/V]')
        grid on;
        box on;
        legend('gain','u(gain)+','u(gain)-');
    end
    
    % generate smooth gain tfer   
    phi(:,1) = cfg.phi_fm.*(f/max(f)).^cfg.phi_pow;
    u_phi(:,1) = (cfg.u_phi_min^2 + (cfg.u_phi_fm.*(f/max(f)).^cfg.phi_pow).^2).^0.5;
    
    if do_plot
        figure;
        semilogx(f,phi)
        hold on;
        semilogx(f,phi + u_phi,'r')
        semilogx(f,phi - u_phi,'r')
        xlabel('f [Hz]')
        ylabel('phase [\mu{}rad]')
        grid on;
        box on;
        legend('phase','u(phase)+','u(phase)-');
    end
    
    % return phase tfer quantity: 
    qu_phi.v = phi;
    qu_phi.u = u_phi;
    
end