clear;
addpath(genpath('../samak'));

runnumbersS = [35410, 35411, 35412, 35413, 35420, 35422];
runnumbersX = [10, 11, 12, 13, 20, 22];
runnumbers = 1:10;


A10 = load(['tritium-bootcamp/runmat',num2str(runnumbersS(1)),'.mat']);
A11 = load(['tritium-bootcamp/runmat',num2str(runnumbersS(2)),'.mat']);
A12 = load(['tritium-bootcamp/runmat',num2str(runnumbersS(3)),'.mat']);
A13 = load(['tritium-bootcamp/runmat',num2str(runnumbersS(4)),'.mat']);
A20 = load(['tritium-bootcamp/runmat',num2str(runnumbersS(5)),'.mat']);
A22 = load(['tritium-bootcamp/runmat',num2str(runnumbersS(6)),'.mat']);



% A1 = load(['tritium-bootcamp/runmat',num2str(runnumbers(1)),'.mat']);
% A2 = load(['tritium-bootcamp/runmat',num2str(runnumbers(2)),'.mat']);
% A3 = load(['tritium-bootcamp/runmat',num2str(runnumbers(3)),'.mat']);
% A4 = load(['tritium-bootcamp/runmat',num2str(runnumbers(4)),'.mat']);
% A5 = load(['tritium-bootcamp/runmat',num2str(runnumbers(5)),'.mat']);
% A6 = load(['tritium-bootcamp/runmat',num2str(runnumbers(6)),'.mat']);
% A7 = load(['tritium-bootcamp/runmat',num2str(runnumbers(7)),'.mat']);
% A8 = load(['tritium-bootcamp/runmat',num2str(runnumbers(8)),'.mat']);
% A9 = load(['tritium-bootcamp/runmat',num2str(runnumbers(9)),'.mat']);
% A10C = load(['tritium-bootcamp/runmat',num2str(runnumbers(10)),'.mat']);

for ii = 1:6
    
    plot(A10.td)
    plot(A11.td)
    plot(A12.td)
    plot(A13.td)
    plot(A20.td)
    plot(A22.td)
end



td00 = A10.td+A11.td+A12.td+A13.td+A20.td+A22.td;
%td0 = A1.td+A2.td+A3.td+A4.td+A5.td+A6.td+A7.td+A8.td+A9.td+A10C.td;

hold on
for ii = 1:6
    plot(A10.td)
    plot(A11.td)
    plot(A12.td)
    plot(A13.td)
    plot(A20.td)
    plot(A22.td)
end
hold off

TimeSec00 = sum(td00);
%TimeSec0 = sum(td0);

qUfrac00 = td00/TimeSec00;
TimeYear00 = TimeSec00/(365.25*24*60*60);

%qUfrac0 = td0/TimeSec0;
%TimeYear0 = TimeSec0/(365.25*24*60*60);

[nqU,npixel] = size(A11.TBDISallPixels);
[nqUC,npixelC] = size(A1.TBDISallPixels);


qu10 = A10.TBDISallPixels(:,1:3:148*3);
qu11 = A11.TBDISallPixels(:,1:3:148*3);
qu12 = A12.TBDISallPixels(:,1:3:148*3);
qu13 = A13.TBDISallPixels(:,1:3:148*3);
qu20 = A20.TBDISallPixels(:,1:3:148*3);
qu22 = A22.TBDISallPixels(:,1:3:148*3);

% qu1 = A1.TBDISallPixels(:,1:3:148*3);
% qu2 = A2.TBDISallPixels(:,1:3:148*3);
% qu3 = A3.TBDISallPixels(:,1:3:148*3);
% qu4 = A4.TBDISallPixels(:,1:3:148*3);
% qu5 = A5.TBDISallPixels(:,1:3:148*3);
% qu6 = A6.TBDISallPixels(:,1:3:148*3);
% qu7 = A7.TBDISallPixels(:,1:3:148*3);
% qu8 = A8.TBDISallPixels(:,1:3:148*3);
% qu9 = A9.TBDISallPixels(:,1:3:148*3);
% qu10C = A10C.TBDISallPixels(:,1:3:148*3);
%% plotting qU distribution 1
qu00 = (qu10+qu11+qu12+qu13+qu20+qu22)/6;
%qu0 = (qu1+qu2+qu3+qu4+qu5+qu6+qu7+qu8+qu9+qu10C)/10;
format long;
pixel = 33;
bin1 = 1; bin2 = 12;
bin3 = 23; bin4 = 36;
quhist1 = [qu10(bin1, pixel), qu11(bin1, pixel), qu12(bin1, pixel),...
    qu13(bin1, pixel), qu20(bin1, pixel), qu22(bin1, pixel)] ;
quhist2 = [qu10(bin2, pixel), qu11(bin2, pixel), qu12(bin2, pixel),...
    qu13(bin2, pixel), qu20(bin2, pixel), qu22(bin2, pixel)] ;
quhist3 = [qu10(bin3, pixel), qu11(bin3, pixel), qu12(bin3, pixel),...
    qu13(bin3, pixel), qu20(bin3, pixel), qu22(bin3, pixel)] ;
quhist4 = [qu10(bin4, pixel), qu11(bin4, pixel), qu12(bin4, pixel),...
    qu13(bin4, pixel), qu20(bin4, pixel), qu22(bin4, pixel)] ;
quhist1 = quhist1-qu00(bin1, pixel);
quhist2 = quhist3-qu00(bin2, pixel);
quhist3 = quhist3-qu00(bin3, pixel);
quhist4 = quhist4-qu00(bin4, pixel);

figure(11);
subplot(2,2,1);
nhist(quhist1);
title(sprintf('Pixel %u Bin %u',pixel, bin1)); 
xlabel('\Delta qU [eV]');
ylabel('');

subplot(2,2,2);
nhist(quhist2);
title(sprintf('Pixel %u Bin %u',pixel, bin2)); 
xlabel('\Delta qU [eV]');
ylabel('');

subplot(2,2,3);
nhist(quhist3);
title(sprintf('Pixel %u Bin %u',pixel, bin3)); 
xlabel('\Delta qU [eV]');
ylabel('');

subplot(2,2,4);
nhist(quhist4);
title(sprintf('Pixel %u Bin %u',pixel, bin4)); 
xlabel('\Delta qU [eV]');
ylabel('');
%% plotting 2 

bin1 = 1; bin2 = 12;
bin3 = 23; bin4 = 36;
quhist10 = [qu10(bin1,:)-qu00(bin1, :), qu11(bin1, :)-qu00(bin1, :), qu12(bin1, :)-qu00(bin1, :),...
    qu13(bin1, :)-qu00(bin1, :), qu20(bin1, :)-qu00(bin1, :), qu22(bin1, :)-qu00(bin1, :)] ;
quhist20 = [qu10(bin2,:)-qu00(bin2, :), qu11(bin2, :)-qu00(bin2, :), qu12(bin2, :)-qu00(bin2, :),...
    qu13(bin2, :)-qu00(bin2, :), qu20(bin2, :)-qu00(bin2, :), qu22(bin2, :)-qu00(bin2, :)] ;
quhist30 = [qu10(bin3, :)-qu00(bin3, :), qu11(bin3, :)-qu00(bin3, :), qu12(bin3, :)-qu00(bin3, :),...
    qu13(bin3, :)-qu00(bin3, :), qu20(bin3, :)-qu00(bin3, :), qu22(bin3, :)-qu00(bin3, :)] ;
quhist40 = [qu10(bin4, :)-qu00(bin4, :), qu11(bin4, :)-qu00(bin4, :),qu12(bin4, :)-qu00(bin4, :),...
    qu13(bin4, :)-qu00(bin4, :), qu20(bin4, :)-qu00(bin4, :), qu22(bin4, :)-qu00(bin4, :)] ;

figure(22);
subplot(2,2,1);
nhist(quhist10);
title(sprintf('All Pixel Bin %u', bin1)); 
xlabel('\Delta qU [eV]');
ylabel('');

subplot(2,2,2);
nhist(quhist20);
title(sprintf('Bin %u', bin2)); 
xlabel('\Delta qU [eV]');
ylabel('');

subplot(2,2,3);
nhist(quhist30);
title(sprintf('Bin %u', bin3)); 
xlabel('\Delta qU [eV]');
ylabel('');

subplot(2,2,4);
nhist(quhist40);
title(sprintf('Bin %u', bin4)); 
xlabel('\Delta qU [eV]');
ylabel('');


%% 


c10 = A10.TBDISallPixels(:,2:3:148*3);
c11 = A11.TBDISallPixels(:,2:3:148*3);
c12 = A12.TBDISallPixels(:,2:3:148*3);
c13 = A13.TBDISallPixels(:,2:3:148*3);
c20 = A20.TBDISallPixels(:,2:3:148*3);
c22 = A22.TBDISallPixels(:,2:3:148*3);

% c1 = A1.TBDISallPixels(:,2:3:148*3);
% c2 = A3.TBDISallPixels(:,2:3:148*3);
% c3 = A3.TBDISallPixels(:,2:3:148*3);
% c4 = A4.TBDISallPixels(:,2:3:148*3);
% c5 = A5.TBDISallPixels(:,2:3:148*3);
% c6 = A6.TBDISallPixels(:,2:3:148*3);
% c7 = A7.TBDISallPixels(:,2:3:148*3);
% c8 = A8.TBDISallPixels(:,2:3:148*3);
% c9 = A9.TBDISallPixels(:,2:3:148*3);
% c10C = A10C.TBDISallPixels(:,2:3:148*3);


c00 = (c10+c11+c12+c13+c20+c22);
%c0 = (c1+c2+c3+c4+c5+c6+c7+c8+c9+c10C);


TBDIS35000 = zeros(nqU,npixel);

TBDIS35000(:,1:3:148*3) = qu00;
TBDIS35000(:,2:3:148*3) = c00;
TBDIS35000(:,3:3:148*3) = sqrt(c00);

TBDIS0 = zeros(nqUC,npixelC);


%TBDIS0(:,1:3:148*3) = qu0;
%TBDIS0(:,2:3:148*3) = c0;
%TBDIS0(:,3:3:148*3) = sqrt(c0);

A00 = A10;

A00.td = td00;
A00.TBDISallPixels = TBDIS35000;
A00.TimeSec = TimeSec00;
A00.qUfrac = qUfrac00;
A00.TimeYear = TimeYear00;

%A0 = A1;

%A0.td = td0;
% A0.TBDISallPixels = TBDIS0;
% A0.TimeSec = TimeSec0;
% A0.qUfrac = qUfrac0;
% A0.TimeYear = TimeYear0;


%save('runmat35410_35422.mat','-struct','A00');

%save('runmat1_10.mat','-struct','A0');





%sum(1,2,3);
