M1=importdata('Measurement1.txt');
M2=importdata('Measurement2.txt');
B =importdata('Bkg.txt');
energy=0.353989.*M1(:,1)-0.037;
energyB=0.35395.*M1(:,1)-0.037;
timeshort=8128.204;
timeBkg=84381.277;
timelong=64433.033;
FilterTime=4182;

for i=1:numel(M1(:,1))-2
    B(i,2) = B(round(0.353989/0.3539461*i),2);
end
M1MinusB = M1(:,2)-B(:,2)*timeshort/timeBkg;
M2MinusB = M2(:,2)-B(:,2)*timelong/timeBkg;
M1MinusB((end-2):end)=0;
M2MinusB((end-2):end)=0;

Ngamma_214Pb_low        = sum(M1(find(abs(energy-294.2)<0.177):find(abs(energy-296.2)<0.177),2)-B(find(abs(energy-294.2)<0.177):find(abs(energy-296.2)<0.177),2).*timeshort./timeBkg);
Ngamma_214Pb_high_short = sum(M1(find(abs(energy-350.9)<0.177):find(abs(energy-352.9)<0.177),2)-B(find(abs(energy-350.9)<0.177):find(abs(energy-352.9)<0.177),2).*timeshort./timeBkg);
Ngamma_214Pb_high_long  = sum(M2(find(abs(energy-350.9)<0.177):find(abs(energy-352.9)<0.177),2)-B(find(abs(energy-350.9)<0.177):find(abs(energy-352.9)<0.177),2).*timelong./timeBkg);
Ngamma_212Pb            = sum(M1(find(abs(energy-237.6)<0.177):find(abs(energy-239)<0.177),2)-B(find(abs(energy-237.6)<0.177):find(abs(energy-239)<0.177),2).*timeshort./timeBkg);

[RhoPb_214_1, RhoRn_222_1]=Densities(Ngamma_214Pb_high_short,efficiency(351.9),0.371,1/183,1/1608,1/330912,0,timeshort,FilterTime);
[RhoPb_214_2, RhoRn_222_2]=Densities(Ngamma_214Pb_high_long,efficiency(351.9),0.371,1/183,1/1608,1/330912,0,timelong,FilterTime);
[RhoPb_214_3, RhoRn_222_3]=Densities(Ngamma_214Pb_low,efficiency(295.2),0.192,1/183,1/1608,1/330912,0,timeshort,FilterTime);
[RhoPb_212, RhoRn_220]=Densities(Ngamma_212Pb,efficiency(238.6),0.436,1/0.15,1/38304,1/55.6,0,timeshort,FilterTime);

plot(energy,M2(:,2));
hold on;
plot(energy,B(:,2)*timelong/timeBkg);
plot(energy,M2MinusB);
hold off;

function DetEff=efficiency(energy)
    DetEff=(1.74123*exp(-0.00553445*energy+1.68653)+0.51814)/100;
end

function [RhoPb, RhoRn]=Densities(Ngamma,epsd,R,lambdaPo,lambdaPb,lambdaRn,t1,t2,Tf)
    V=1/18; %m³/s
    epsf=0.9997;
    RhoPb=Ngamma/(epsf*epsd*R*V*((1/lambdaPo+1/lambdaPb+lambdaPb/(lambdaPo*(lambdaPo-lambdaPb))+lambdaPo*exp(-lambdaPb*Tf)/(lambdaPb*(lambdaPb-lambdaPo)))*...
        (exp(-lambdaPb*t1)-exp(-lambdaPb*t2))+...
        (lambdaPb/(lambdaPo*(lambdaPb-lambdaPo))*(1-exp(-lambdaPo*Tf)))*(exp(-lambdaPo*t1)-exp(-lambdaPo*t2))));
    RhoRn=lambdaPb*RhoPb/lambdaRn;
end