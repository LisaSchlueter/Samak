function [ pvalue ] = chi2pvalue( x, ndof )
% chi2cdf with a different name
% Input: 
%  - x    : chi^2
%  - ndof : number of degrees of freedom
%
% Output:
%  - p-value
%

if ndof < 0
    pvalue=0
end

if x<=0
    if x<0
        pvalue=0
    else
        pvalue=1
    end
end

pvalue=gammainc(x/2,ndof/2,'upper');
    
end

