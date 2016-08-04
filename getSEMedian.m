% [se,bs] = getSEMedian(X,N) returns the standard error of the median of
% the values in X. 

function [se,bs] = getSEMedian(X,N)

if ~exist('N','var');                   N=length(X);                    end

bs = bootstrp(N,@median,X);
se = std(bs);
end