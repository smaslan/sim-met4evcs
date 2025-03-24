function [num_str,num_with_unit,fmt,si_prefix,si_scale] = num_fmt(num, abs_unc, rel_unc, num_width, si_max, unit)
% Simple number formater by uncertainty/accuracy. Can do automatic SI prefix scaling.
%
% Usage:
%   [num_str,num_with_unit,fmt,si_prefix,si_scale] = num_fmt(num, abs_unc, rel_unc, num_width)
%   [num_str,num_with_unit,fmt,si_prefix,si_scale] = num_fmt(num, abs_unc, rel_unc, num_width, si_max)
%   [num_str,num_with_unit,fmt,si_prefix,si_scale] = num_fmt(num, abs_unc, rel_unc, num_width, si_max, unit)
%
% Example:
%   [num_str,num_with_unit,fmt,si_prefix,si_scale] = num_fmt(12.34567e3, 0.0001, 100e-6, 8, 100, 'Wh')
%   num_str = ' 12.3457'
%   num_with_unit = ' 12.3457 kWh'
%   fmt = '%8.4f';
%   si_prefix = 'k'
%   si_scale = 0.001
%
% Parameters:
%   num - value to print
%   abs_unc - absolute uncertainty
%   rel_unc - relative uncertainty [-]
%   num_width - print characters count for entire number (including decimal dot)
%   si_max - optional, auto SI prefix scaling rule (by default disabled), e.g.:
%            si_max = 100 for num = 99 will print '99 [unit]',
%            for num = 101 will print '0.101 k[unit]'
%   unit - optional, quantity unit string
%
% Returns:
%   num_str - formatted number string (no unit or SI prefix)
%   num_with_unit - formatted number string with SI prefix and unit
%   fmt - number format string
%   si_prefix - SI prefix character or '' for none
% 
% This is part of the EVCS charging waveform simulator.
% Developed in scope of EPM project 23IND06 Met4EVCS: https://www.vsl.nl/en/met4evcs/
% (c) 2024, Stanislav Maslan (smaslan@cmi.cz)
% The script is distributed under MIT license, https://opensource.org/licenses/MIT.
%

    sip = 5;
    tmp = num; 
    if exist('si_max','var')
        while abs(tmp) && abs(tmp) > si_max
            tmp = tmp*0.001;
            sip = sip + 1;
        end
        while abs(tmp) && abs(tmp) < 0.001*si_max
            tmp = tmp*1000.0;
            sip = sip - 1;            
        end
    end
    
    sip_list = 'pnum kMGTP';
    sip = min(max(sip,1),numel(sip_list));
    si_scale = 10^((5 - sip)*3); 
    num = num*si_scale;
    si_prefix = sip_list(sip);
    if si_prefix == ' '
        si_prefix = '';
    end 
    
    acc = max(abs(num)*rel_unc, abs_unc*si_scale);
    
    frac = max(ceil(-log10(acc))+1,0);
    
    fmt = sprintf('%%%.0f.%.0ff',num_width,frac);
        
    if ~exist('unit','var')
        unit = '';
    end
    unit = [si_prefix unit];
    
    num_str = sprintf([fmt], num);
    num_with_unit = sprintf([fmt ' %s'], num, unit);

end