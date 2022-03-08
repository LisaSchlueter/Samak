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

% [linFitpar3, linFiterr3, ~,~] = linFit(R.RingList(1:3)',mNuSq(1:3),mNuSqErr(1:3));
% linFitpar3(1)./linFiterr3(1)
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

%% fit mNuSq as a function of actual radial FPD position
% radial ring position (from FPD Viewer)
rStart = [0 0.7398 1.4796 1.9573 2.3394 2.6674 2.9592 3.2247 ...
    3.4699 3.699 3.9146 4.119 4.3137]';
rEnd = [0.7398 1.4796 1.9573 2.3394 2.6674 2.9592 3.2247 ...
    3.4699 3.699 3.9146 4.119 4.3137 4.5]';
Radius_cm = mean([rStart';rEnd']);
if strcmp(A.RingMerge,'Full')
Radius_Ring_cm = zeros(4,1); %average radial position of pseudo-ring
Radius_Ring_cm(1) = mean(Radius_cm(1:3));
Radius_Ring_cm(2) = mean(Radius_cm(4:6));
Radius_Ring_cm(3) = mean(Radius_cm(7:9));
Radius_Ring_cm(4) = mean(Radius_cm(10:12));
end
[linFitpar, linFiterr, linFitchi2min,linFitdof] = linFit(Radius_Ring_cm,mNuSq,mNuSqErr);
fprintf('Linear fit mNuSq(radius): slope = (%.1f +- %.1f) eV2/cm -> %.1f sigma \n',linFitpar(1),linFiterr(1),abs(linFitpar(1)./linFiterr(1)));





