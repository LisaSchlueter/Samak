% uniform fit on knm2 stacked data
% settings
RunList = 'KNM2_RW1';
fixPar = 'E0 Bkg Norm';% free parameter%'1 5 6 7 8 9 10 11 12'; % fixed parameter
DataType = 'Real';
FSDFlag = 'Sibille0p5eV';
NonPoissonScaleFactor = 1.064;
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
exclDataStart = 2; % 12 == 40 eV range
chi2 = 'chi2Stat';
RunAnaArg = {'RunList',RunList,'fixPar',fixPar,'DataType',DataType,...
            'FSDFlag',FSDFlag,'ELossFlag',ELossFlag,...
            'NonPoissonScaleFactor',NonPoissonScaleFactor,'exclDataStart',exclDataStart,...
            'AnaFlag',AnaFlag,'chi2',chi2};

% read data and set up model
A = MultiRunAnalysis(RunAnaArg{:});

A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
TBDDS_i = A.ModelObj.TBDDS;
TBDIS_i = A.ModelObj.TBDIS./(A.ModelObj.qUfrac.*A.ModelObj.TimeSec);
Te_i = A.ModelObj.Te;
%%
A.ModelObj.DopplerEffectFlag = 'matConv';
A.ModelObj.ComputeTBDDS;A.ModelObj.ComputeTBDIS;
TBDDS_Doppler = A.ModelObj.TBDDS;
TBDIS_Doppler = A.ModelObj.TBDIS./(A.ModelObj.qUfrac.*A.ModelObj.TimeSec);
Te_Doppler = A.ModelObj.Te;
%% differential spectra
pDS1 = plot(A.ModelObj.Te-A.ModelObj.Q,TBDDS_i,'LineWidth',3);
hold on;
pDS2 = plot(A.ModelObj.Te-A.ModelObj.Q,TBDDS_Doppler,'LineWidth',pDS1.LineWidth);
PrettyFigureFormat;
hold off;
leg = legend([pDS1,pDS2],'Without Doppler effect','With Doppler effet');
legend boxoff
leg.Location='southwest';
xlabel(sprintf('{\\itE} - {\\itE}_0 (eV)'));
ylabel(sprintf('d\\Gamma/dE (cps)'))
xlim([-15,0]);
set(gca,'YScale','log');

%% integral spectra: rel. difference
pIS1 = plot(A.ModelObj.qU(2:end)-18575,100.*(TBDIS_Doppler(2:end)-TBDIS_i(2:end))./TBDIS_i(2:end),'LineWidth',3);
PrettyFigureFormat('FontSize',24);
leg = legend('(Doppler - No Doppler)/ No Doppler ');
legend boxoff
%leg.Location='southwest';
xlabel(sprintf('{\\itE} - {\\itE}_0 (eV)'));
ylabel(sprintf('Rel. difference (%%)'))
xlim([-15,2]);
set(gca,'YScale','lin');

% plot(A.ModelObj.qU-18575,(TBDIS_i));
% hold on
% plot(A.ModelObj.qU-18575,(TBDIS_Doppler));