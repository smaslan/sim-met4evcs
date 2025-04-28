function [f,qu_gain,qu_phi] = gen_shunt_tfer(f_max,f_count, cfg, varargin)
% Simple function to simulate parametric frequency transfer of current shunt.
% Uses model: Z = (Rs + Ls) || Cp
%
% Usage:
%   [f,qu_gain,qu_phi] = gen_shunt_tfer(f_max,f_count, cfg) 
%   [f,qu_gain,qu_phi] = gen_shunt_tfer(f_max,f_count, cfg, debug_plot)
%
% Parameters:
%   f_max - max frequency of freq transfer (usually nyquist fs/2)
%   f_count - number of frequency spots to generate (logspace)
%   cfg - transfer parametrs:
%   cfg.Rs - series resistance [Ohm]
%   cfg.u_Rs - absolute Rs uncertainty [Ohm]
%   cfg.Ls - series inductance [H]
%   cfg.u_Ls - absolute Ls uncertainty [H]
%   cfg.Cp - parallel capacitance [F]
%   cfg.u_Cp - absolute Cp uncertainty [F]
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
% Source: https://github.com/smaslan/sim-met4evcs
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.

    f = logspaced(10,f_max+eps,f_count);
    f = [0;f(:)];
    f(end) = f_max + eps;
    w = 2*pi*f;
    
    mcc = 10000;
    Rs = cfg.Rs + randn(1,mcc)*cfg.u_Rs;
    Ls = cfg.Ls + randn(1,mcc)*cfg.u_Ls;
    Cp = cfg.Cp + randn(1,mcc)*cfg.u_Cp;
    
    Z = 1./(1./(Rs + j*w.*Ls) + j*w.*Cp);
    
    gain = abs(Z);
    phi = angle(Z);
    
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
        ylabel('|Z| [\Omega]')
        grid on;
        box on;
        legend('|Z|','u(|Z|)+','u(|Z|)-');
        
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