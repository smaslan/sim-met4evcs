function [fx_rel, hx_rel, phx] = gen_thd_harms(thd, shape, max_n)
% Generate harmonic components to match desired THD value.
%
% Usage:
%  [fx_rel, hx_rel, phx] = gen_thd_harms(thd, shape, max_n)
%
% Parameters:
%  thd - desired THD value [%]
%  shape - shape of harmonic spectrum:
%          'exp' - exponentially decaying harmonics
%          'sqr' - square wave harmonics (odd harmonics)
%  max_n - maximum harmonic index to generate, e.g. 5 will generate harmonics up to 5th
%
% Returns:
%  fx_rel - relative frequencies (indices) of harmonics, e.g. [3, 5, 7]
%  hx_rel - harmonic amplitudes relative to fundamental, e.g. [0.5, 0.3, 0.2]
%  phx - phase angles of harmonics [rad]
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% Source: https://github.com/smaslan/sim-met4evcs
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.    

    if strcmpi(shape, 'exp')
        
        fx_rel = [2:max_n];
        hx = exp(1./fx_rel)/exp(1/fx_rel(1))*0.01*thd;
        thd_x = sum(hx.^2).^0.5;    
        hx_rel = hx*0.01*thd/thd_x;
        phx = repmat(0, size(hx_rel));
        
    elseif strcmpi(shape, 'sqr')
        
        fx_rel = [3:2:max_n];
        hx = 1./fx_rel*0.01*thd;
        thd_x = sum(hx.^2).^0.5;    
        hx_rel = hx*0.01*thd/thd_x;
        hx_rel = hx;
        phx = repmat(0, size(hx_rel));
    
    else
        error(sprintf('Unknown THD shape ''%s''!', shape));
    end

end