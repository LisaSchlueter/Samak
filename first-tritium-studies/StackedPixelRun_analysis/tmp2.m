x = 1:100;
y = x.^2;
yd = gradient(y);

figure(999)
plot(x,y)
hold on
plot(x,yd)
hold off