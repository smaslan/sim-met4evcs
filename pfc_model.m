function [f_pfc] = pfc_model(phi, f_min, f_max)
% very crude model of crm pfc switching frequency vs voltage phase phi
% f_min and f_max are minimum and maximum switching frequnecies     
    f = (mod(phi,2*pi) - pi).^2;
    f_pfc = f/max(f)*(f_max - f_min) + f_min;    
end


