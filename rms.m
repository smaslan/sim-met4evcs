function [rms] = rms(A)
% get RMS value of vector    
    rms = sum(A.^2).^0.5;
end