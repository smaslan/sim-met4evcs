function [x] = rand2(sz,mid,rng)
% generate random number(s) x = mid +- rng
    x = (rand(sz)*2 - 1).*rng + mid;	
end