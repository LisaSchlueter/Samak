% uniform fit on knm2 stacked data
% settings
RunList = 'KNM1';
fixPar = 'mNu E0 Bkg Norm'; % free parameter
DataType = 'Real';
FSDFlag = 'Sibille0p5eV';
NonPoissonScaleFactor = 1.064;
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
range = 40; 
chi2 = 'chi2Stat';
RingMerge = 'Full';%None';
RunAnaArg = {'RunList',RunList,...
             'fixPar',fixPar,...
             'DataType',DataType,...
            'FSDFlag',FSDFlag,...
            'ELossFlag',ELossFlag,...
            'NonPoissonScaleFactor',NonPoissonScaleFactor,...
            'AnaFlag',AnaFlag,...
            'chi2',chi2,...
            'RingMerge',RingMerge,...
            'AngularTFFlag','OFF'};

% read data and set up model
A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart = A.GetexclDataStart(range);
R = RingAnalysis('RunAnaObj',A,'RingList',A.RingList);

%%
R.FitRings('SaveResult','ON','RecomputeFlag','OFF');

%%
close all
R.PlotFits('SavePlot','ON','Blind','OFF','PlotPar',1,'PlotMode','Abs','YLim',[-8.5 6.5])
%R.PlotFits('SavePlot','ON','Blind','OFF','PlotPar',2,'PlotMode','Abs','YLim',18573+[+0.25,+1.15])
%R.PlotFits('SavePlot','ON','Blind','OFF','PlotPar',3,'PlotMode','Abs','YLim',[2.1 3.1]);%,'YLim',[-8.5 6.5])
%R.PlotFits('SavePlot','ON','Blind','OFF','PlotPar',4,'PlotMode','Rel','YLim',[-1.3e-2 +2.2e-2]);%,'YLim',[-8.5 6.5])%


%%
mNuSq = R.FitResult.par(:,1);
mNuSqErr = 0.5*(abs(R.FitResult.errNeg(:,1))+R.FitResult.errPos(:,1));
               
R.FitResult.err(:,1);

% compare to constant hypothesis
mean = wmean(mNuSq,1./mNuSqErr.^2);
chi2mean = sum((mNuSq-mean).^2./mNuSqErr.^2);

%
% compare to linear hypothesis
[linFitpar, linFiterr, linFitchi2min,linFitdof] = linFit(R.RingList',mNuSq,mNuSqErr);
linFitpar(1)./linFiterr(1)

[linFitpar3, linFiterr3, ~,~] = linFit(R.RingList(1:3)',mNuSq(1:3),mNuSqErr(1:3));
linFitpar3(1)./linFiterr3(1)
%%
Deltachi2 = chi2mean-linFitchi2min;
chi2cdf(Deltachi2,1)

%1.20
%% some statistics
BkgRate_r = arrayfun(@(x) x.ModelObj.BKG_RateSec_i,R.MultiObj)+R.FitResult.par(:,4);
Bkg_r = BkgRate_r.*A.ModelObj.qUfrac(A.exclDataStart:end)'.*A.ModelObj.TimeSec;
TBDIS_r = cell2mat(arrayfun(@(x) x.RunData.TBDIS(A.exclDataStart:end),R.MultiObj,'UniformOutput',false)')';
Signal_r = TBDIS_r-Bkg_r;
sum(Signal_r')



