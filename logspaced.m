function [y] = logspaced(a,b,n,digits)
% Generate logspace() from a to b, wit n steps and optionally round to digits
	
    y = logspace(log10(a),log10(b),n);
    
    if exist('digits','var') && digits          
        % round to digits count
        digz = ceil(log10(y));    
        round_base = 10.^-(digz - digits);    
        y = round(y.*round_base)./round_base;
    end
    y = unique(y);
    
end