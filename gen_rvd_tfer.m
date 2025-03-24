function [f,qu_gain,qu_phi] = gen_rvd_tfer(f_max,f_count, cfg, varargin)
% Simple function to simulate parametric frequency transfer of resistive voltage divider.
%
% Usage:
%   [f,qu_gain,qu_phi] = gen_rvd_tfer(f_max,f_count, cfg) 
%   [f,qu_gain,qu_phi] = gen_rvd_tfer(f_max,f_count, cfg, debug_plot)
%
% Parameters:
%   f_max - max frequency of freq transfer (usually nyquist fs/2)
%   f_count - number of frequency spots to generate (logspace)
%   cfg - transfer parametrs:
%   cfg.Rhi - high-side parallel resistance [Ohm]
%   cfg.u_Rhi - absolute Rhi uncertainty [Ohm]
%   cfg.Chi - high-side parallel capacitance [F]
%   cfg.u_Chi - absolute Chi uncertainty [F]
%   cfg.Rlo - low-side parallel resistance [Ohm]
%   cfg.u_Rlo - absolute Rlo uncertainty [Ohm]
%   cfg.Clo - low-side parallel capacitance [F]
%   cfg.u_Clo - absolute Clo uncertainty [F]
%   cfg.Llo - low-side series inductance [H]
%   cfg.u_Llo - absolute Llo uncertainty [H]
%   debug_plot - optional plot of generated tfers
%
% Returns:
%   f - generated frequency vector from 0 to f_max [Hz]
%   qu_gain.v - gain transfer (modulus of impedance) [Ohm]   
%   qu_gain.u - absolute gain transfer uncertainty [Ohm]
%   qu_phi.v - phase transfer [rad]   
%   qu_phi.u - absolute phase transfer uncertainty [rad]
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.

    f = logspaced(10,f_max+eps,f_count);
    f = [0;f(:)];
    f(end) = f_max + eps;
    w = 2*pi*f;
    
    mcc = 10000;
    Rhi = cfg.Rhi + randn(1,mcc)*cfg.u_Rhi;
    Chi = cfg.Chi + randn(1,mcc)*cfg.u_Chi;
    Rlo = cfg.Rlo + randn(1,mcc)*cfg.u_Rlo;
    Clo = cfg.Clo + randn(1,mcc)*cfg.u_Clo;
    Llo = cfg.Llo + randn(1,mcc)*cfg.u_Llo;
    
    Zhi = 1./(1./Rhi + j*w.*Chi);
    Zlo = 1./(1./Rlo + j*w.*Clo) + j*w.*Llo;
       
    tf = Zlo./(Zlo + Zhi);
    
    gain = abs(tf);
    phi = angle(tf);
    
    qu_gain.v = mean(gain,2);
    qu_gain.u = std(gain,[],2);
    qu_phi.v = mean(phi,2);
    qu_phi.u = std(phi,[],2);
    
    % debug plot?
    do_plot = nargin() > 3 && isnumeric(varargin{1}) && varargin{1};
    
    if do_plot
        figure;
        semilogx(0.001*f(2:end),qu_gain.v(2:end))
        hold on;
        semilogx(0.001*f(2:end),qu_gain.v(2:end) + qu_gain.u(2:end),'r')
        semilogx(0.001*f(2:end),qu_gain.v(2:end) - qu_gain.u(2:end),'r')
        xlabel('f [kHz]')
        ylabel('gain [V/V]')
        grid on;
        box on;
        legend('gain','u(gain)+','u(gain)-');
        
        figure;
        semilogx(0.001*f(2:end),qu_phi.v(2:end))
        hold on;
        semilogx(0.001*f(2:end),qu_phi.v(2:end) + qu_phi.u(2:end),'r')
        semilogx(0.001*f(2:end),qu_phi.v(2:end) - qu_phi.u(2:end),'r')
        xlabel('f [kHz]')
        ylabel('\Phi [rad]')
        grid on;
        box on;
        legend('\Phi','u(\Phi)+','u(\Phi)-');
        
        
        tau = qu_phi.v./w;
        u_tau = qu_phi.u./w;
        
        figure;
        semilogx(0.001*f(2:end),1e9*tau(2:end))
        hold on;
        semilogx(0.001*f(2:end),1e9*(tau(2:end) + u_tau(2:end)),'r')
        semilogx(0.001*f(2:end),1e9*(tau(2:end) - u_tau(2:end)),'r')
        xlabel('f [kHz]')
        ylabel('\tau [ns]')
        grid on;
        box on;
        legend('\tau','u(\tau)+','u(\tau)-');
        
    end 
    

end