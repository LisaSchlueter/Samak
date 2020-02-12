
addpath(genpath('../../../Samak2.0'));
PixelID = 1;
mnuSq_t = 0;
S0 = init_bootcampR35410('PixelID',PixelID,'mnuSq_i',mnuSq_t);
D0 = load('data/runmat35410.mat');
S1 = init_bootcampR35411('PixelID',PixelID,'mnuSq_i',mnuSq_t);
D1 = load('data/runmat35411.mat');
S2 = init_bootcampR35412('PixelID',PixelID,'mnuSq_i',mnuSq_t);
D2 = load('data/runmat35412.mat');
S3 = init_bootcampR35413('PixelID',PixelID,'mnuSq_i',mnuSq_t);
D3 = load('data/runmat35413.mat');
S20 = init_bootcampR35420('PixelID',PixelID,'mnuSq_i',mnuSq_t);
D20 = load('data/runmat35420.mat');
S22 = init_bootcampR35422('PixelID',PixelID,'mnuSq_i',mnuSq_t);
D22 = load('data/runmat35422.mat');
quIndex      = (PixelID-1)*3+1;  %[1:3:148*3];
countIndex   = (PixelID-1)*3+2; %[2:3:148*3];
coutErrIndex = (PixelID-1)*3+3; %[3:3:148*3];
Counts0      = D0.TBDISallPixels(:,countIndex);
CountErrors0 = D0.TBDISallPixels(:,coutErrIndex); CountErrors0(CountErrors0==0)=1;
Data0 = [...
    S0.qU,...
    Counts0,...
    CountErrors0];


Run = [35410, 35411 35412 35413 35420 35422];
for i=1:6
[parA, err, chi2min, ndof] =  fit_tbdpixel('Run',Run(i));
ParnPixels(i,:) = [i parA err chi2min ndof];
end
figure(111);
hdata = errorbar(Data0(:,1),Data0(:,2),Data0(:,3),...
    'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
hfit1 = plot(Data0(:,1),tbd_modelint(parA,S0),...
    'Color','Black','LineWidth',1,'LineStyle','-');

hfit3 = line([Data0(1,1) Data0(end,1)],[(S0.BKG_RateSec_i+par(3))*S.TimeSec*S.qUfrac(end) (S.BKG_RateSec_i+par(3))*S.TimeSec*S.qUfrac(end)],'LineStyle','--','Color','Red');
hold off
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Counts/qU','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
%mydata = sprintf('Data: m=%.2f eV - E0=%.2f eV \n',...
%    sqrt(mnuSq_t),18575);
%myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n m=%.3f  +- %.3f eV \n E0=%.3f  +- %.3f meV \n B=%.3f  +- %.3f cps \n N=%.3f  +- %.3f',...
 %   chi2min,S0.nqU-4,sqrt(abs(S.mnuSq_i+par(1))),err(1),(par(2))*1e3,err(2)*1e3,...
 %   (S.BKG_RateSec_i+par(3))*1e3,err(3)*1e3,par(4),err(4));
%mychi2 = sprintf('Data');
%l1 = legend([hdata  hfit1 hfit3],mychi2,myfit,'Background','Location','NorthEast') ;
%l1.FontSize = 11;
axis([min(S0.qU) max(S0.qU)+1 0.*min(Data0(:,2)) max(Data0(:,2))*1.2])
title(sprintf('Pixel %g - Data and Fit',S0.FPD_Pixel),'FontSize',14);