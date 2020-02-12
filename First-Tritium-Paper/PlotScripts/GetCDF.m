function CDF = GetCDF(x,y)
% calculate cumulative distribution functon
CDF = zeros(numel(x),1);
for i=2:numel(x)
    CDF(i) = simpsons(x(1:i),y(1:i));   
end

CDF = CDF./simpsons(x,y);
end