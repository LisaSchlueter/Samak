function out = Intconv(x,xInt,f1,f2)
% numerical convolution of 2 vectors
% alternative to conv.m 
% difference: integration instead of sum

% input:
% 1) x     evaluation vector
% 2) xInt  integration vector
% 3) f1    function handle 1
% 4) f1    function handle 2

% output = convolution evaluted at points x

%% make sure x vectors are in correct format
% xInt==row vector , x== column vector
if isrow(xInt) && isrow(x)
    x = x';
elseif iscolumn(xInt)
    xInt = xInt';
    if isrow(x)
        x = x';
    end
end
%% make sure f are in corrector fomat
% f1 --> column vector
% f2 --> row vector
if isrow(f1(1:2))
    f1 = @(e) f1(e)';
end

if iscolumn(f2(1:2))
    f2 = @(e) f2(e)';
end
%% calc
% f1--> row
Intfun = f1(xInt)'.*f2(x-repmat(xInt,numel(x),1));
out = simpsons(xInt',Intfun')';

end