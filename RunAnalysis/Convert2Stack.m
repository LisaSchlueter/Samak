function y = Convert2Stack(x,startValue)
%x is array to be converted
y = zeros(1,numel(x));
y(1)= startValue;
for i=2:numel(y)
    y(i)= x(i)-sum(y(1:i-1));
end
end