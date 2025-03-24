function [ip] = gen_pfc_emi(t, phi, f_min, f_max, fh_count, fh_max_f, amp)
% Very crude generator of harmonic spurs produced by active PFC circuit in Critical Conduction Mode.
% It generates triangular wave with frequency in range f_min to f_max depending on phase angle phi.
% Maximum harmonics count fh_count and maximum harmonic frequency fh_max_f can be defined.
% Amplitude can be scalar or vector (one element per element of t and phi).
%
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.

    % first sample with positive time
    t0id = find(t >= 0.0,1);     
    dt = abs(t(2) - t(1));
    
    % make amplitude vector if only scalar provided
    if isscalar(amp)
        amp = repmat(amp, size(t));    
    end
    amp(1:t0id-1) = 0;
    
    % generate some pfc model
    f_pfc = pfc_model(phi, f_min, f_max);
    
    % for each harmonic (one by one to save memory):
    ip = zeros(size(t));
    for h = 1:2:fh_count
        
        % this harmonic component frequencies
        fh = f_pfc*h;
        fh(1:t0id-1) = 0;
        
        % this harmonic amplitude
        A = amp/h^2;
        % limit harmonics by frequency
        A(fh > fh_max_f) = 0;
        
        % add one harmonic to mix
        wh = 2*pi*fh;
        ip = ip + A.*cos(cumsum(wh*dt));    
    end
    
end