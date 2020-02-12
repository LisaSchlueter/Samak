function [mnuSq_i_Fit, par, err, chi2min, dof,mNu90,mNumin] = NuMassScan_SensitivityNominal(varargin)
% -----------------------------------------------------------------------------------------------
% Goal: 
% Obtain sensitivity on neutrino mass squared
% Method: Scan
% Simulate KATRIN 1 time (asimov)
% Options: MTD, Background, Time
% Fit with statistic only or with a covariance matrix (E0,N,B free) and a fixed neutrino mass (scan over mnu)
% Plot chi^2 as a function of neutrino mass
% Compute sensitivity to neutrino mass squared
% ------------------------------------------------------------------------------------------------
addpath(genpath('../../../Samak2.0'));
p = inputParser;
p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x) && x>0);
p.addParameter('range',60,@(x)isfloat(x));
p.addParameter('TD','DR30',@(x)ischar(x)); %only for non SensitivityBa....
p.addParameter('BKG_RateSec','',@(x)isfloat(x)); %if empty, filled according to Ba
p.addParameter('FPD_MeanEff',0.9,@(x)isfloat(x) && x>=0);
%Simulation Model Parameter
p.addParameter('MACE_Ba_T',7e-04,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_B_T',3.6*0.7,@(x)isfloat(x) && x>=0);
%Fit Model Paramater
p.addParameter('mNuStart',0.2,@(x)all(isfloat(x))); % Start mNu for neutrino mass scan (eV)
p.addParameter('mNuStop',0.6,@(x)all(isfloat(x))); % Stop mNu for neutrino mass scan (eV)
p.addParameter('ScanPrcsn',0.02,@(x)all(isfloat(x))); % Scan Precision (Acceptable Delta chi2)
p.addParameter('Q_i',18575,@(x)all(isfloat(x))); % To give Bias of Work Function (in sample spectra)
%Fit Settings
p.addParameter('chi2','chi2CM',@(x)ismember(x,{'chi2Stat','chi2CM','chi2CMShape'}));
p.addParameter('SysEffect','FSD'); %1 sys effect e.g. FSD or struct for more than 1
p.addParameter('plotFit','ON',@(x)ismember(x,{'ON','OFF'}));
%Systematics settings
p.addParameter('PlotCM','OFF',@(x)ismember(x,{'ON','OFF'})); %produces a set of plots and a text file with information
p.addParameter('SysBudget','08',@(x)ischar(x)); % defines the covmat reference file
p.addParameter('Scan','ON',@(x)ismember(x,{'ON','OFF'}));   % Off = no scan, 1 asimov fit with free nu-mass
p.addParameter('BkgCM','OFF',@(x)ismember(x,{'ON','OFF'})); % Background Covariance matrix
p.addParameter('Anchor6G','OFF',@(x)ismember(x,{'ON','OFF'})); % Background Option
p.addParameter('Anchor6GValue',335e-3,@(x)isfloat(x)); % 335=FT, [26keV, 32 keV]ROI
p.parse(varargin{:});

TD           = p.Results.TD;
TimeSec      = p.Results.TimeSec;
range        = p.Results.range;
MACE_Ba_T    = p.Results.MACE_Ba_T;
WGTS_B_T     = p.Results.WGTS_B_T;
BKG_RateSec  = p.Results.BKG_RateSec;
Q_i          = p.Results.Q_i;
chi2         = p.Results.chi2;
SysEffect    = p.Results.SysEffect;
plotFit      = p.Results.plotFit;
PlotCM       = p.Results.PlotCM;
SysBudget    = p.Results.SysBudget;
mNuStart     = p.Results.mNuStart;
mNuStop      = p.Results.mNuStop;
ScanPrcsn    = p.Results.ScanPrcsn;
BkgCM        = p.Results.BkgCM;
Scan         = p.Results.Scan;
FPD_MeanEff  = p.Results.FPD_MeanEff; 
Anchor6G     = p.Results.Anchor6G;
Anchor6GValue = p.Results.Anchor6GValue;

if strcmp(chi2,'chi2Stat')
    SysEffect = '';
end

nFitMax = 20;
%---------------------------------------------parser end -----------------------------------------------%
%% Init:
% Model genereate Asimov Data Set (S)
% Model genereate to Fit (M)
% Covariance Matrix fpr Fit
if isempty(TD)
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
end
%Background
if isempty(BKG_RateSec)
BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Anchor6G',Anchor6G,'Anchor6GValue',Anchor6GValue);
end

% Fit settings
fitter = 'minuit';
exclDataStart = 1;

if strcmp(chi2,'chi2Stat')
    SysEffect = '';  
end
%Simulation Model
ModelArg = {'BKG_RateAllFPDSec',BKG_RateSec,'MACE_Ba_T',MACE_Ba_T,...
                           'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',(WGTS_B_T/3.6)*6,...
                           'TimeSec',TimeSec,'TD',TD,'mnuSq_i',0,'Q_i',Q_i,...
                           'FPD_MeanEff',FPD_MeanEff};
S = ref_TBD_NominalKATRIN(ModelArg{:});                    
S.ComputeTBDDS;
S.ComputeTBDIS;
Data = [S.qU, S.TBDIS, sqrt(S.TBDIS)];
%Fit Model
M = ref_TBD_NominalKATRIN(ModelArg{:});
M.ComputeTBDDS; M.ComputeTBDIS;
TBDIS_Fit = zeros(M.nqU,nFitMax); %Save Fit spectra

%Covariance Matrix
if ~strcmp(chi2,'chi2Stat')
    [~, MultiCM, MultiCMFrac, MultiCMShape, MultiCMNorm] = ...
        ref_CovarianceMatrix_NominalKATRIN('RecomputeFlag','OFF','ModelObj',M,'SysEffect',SysEffect,'PlotCM',PlotCM,'SysBudget',SysBudget);
    
    %Init Fit Model again
    M = ref_TBD_NominalKATRIN(ModelArg{:});
    M.ComputeTBDDS; M.ComputeTBDIS;
    if strcmp(BkgCM,'ON')
        [BkgCM, BkgCMFrac, BkgCMShape, BkgCMFracShape] = ComputeCM_Background('StudyObject',M,'plotFit',PlotCM,'nTrials',1000);
        MultiCM      = MultiCM + BkgCM;
        MultiCMShape = MultiCMShape + BkgCMShape;
    end
elseif strcmp(chi2,'chi2Stat')
    if strcmp(BkgCM,'OFF')
        MultiCM      = 0;
        MultiCMFrac  = 0;
        MultiCMShape = 0;
        MultiCMNorm  = 0;
    elseif strcmp(BkgCM,'ON')
        [BkgCM, BkgCMFrac, BkgCMShape, BkgCMFracShape] = ComputeCM_Background('StudyObject',M,'plotFit',PlotCM,'nTrials',1000);
        MultiCM      = BkgCM;
        MultiCMFrac  = 0;
        MultiCMShape = BkgCMShape;
        MultiCMNorm  = 0;
    end
end
%% Scan
if strcmp(Scan,'ON')
tic
mnuSq_i_Fit          = zeros(nFitMax,1); % neutrino mass squared scan values
mnuSq_i_Fit(2)       = mNuStart^2;    % start value for scan
mnuSq_i_Fit(nFitMax) = mNuStop^2;   % stop value for scan
F                    = cell(length(nFitMax),1);
par                  = zeros(6,nFitMax);
err                  = zeros(6,nFitMax);
chi2min              = zeros(nFitMax,1);
dof                  = zeros(nFitMax,1);

% Test of mNu=0 and Stop Value
% mNu=0 should yield chi2=0
% Stopvalue should yield chi2>2.7
TestVal = [1,nFitMax];
for i=1:numel(TestVal)
    M.mnuSq_i = mnuSq_i_Fit(TestVal(i)); %zero
    F{TestVal(i)} = FITC('SO',M,'DATA',Data,'fitter',fitter,...
        'chi2name',chi2,...
        'COVMAT', MultiCM+diag(S.TBDIS),'COVMATFrac', MultiCMFrac+diag(1./S.TBDIS),...
        'COVMATShape', MultiCMShape,'COVMATNorm',MultiCMNorm,...
        'fixPar','1 5 6',...
        'exclDataStart',exclDataStart);
    par(:,TestVal(i)) = F{TestVal(i)}.RESULTS{1};
    par(1,TestVal(i)) = par(1,TestVal(i))+M.mnuSq_i;
    err(:,TestVal(i)) = F{TestVal(i)}.RESULTS{2};
    chi2min(TestVal(i)) = F{TestVal(i)}.RESULTS{3};
    dof(TestVal(i)) = F{TestVal(i)}.RESULTS{5};
    TBDIS_Fit(:,TestVal(i)) = F{TestVal(i)}.SO.TBDIS;
    if  chi2min(1)>1e-05 %if chi2 isnt 0, problem!
        fprintf(2,'Error: \chi^2 is not 0 for neutrino mass 0! \n');
        return
    elseif TestVal(i)==nFitMax && chi2min(nFitMax)<2.7
        fprintf(2,'Error: Neutrino Mass Scan Range is too small. Enlarge mNuStop! \n');
        return
    end
end

% Do Scan
for i=2:nFitMax
    M.mnuSq_i = mnuSq_i_Fit(i);
    F{i} = FITC('SO',M,'DATA',Data,'fitter',fitter,...
        'chi2name',chi2,...
        'COVMAT', MultiCM+diag(S.TBDIS),'COVMATFrac', MultiCMFrac+diag(1./S.TBDIS),...
        'COVMATShape', MultiCMShape+diag(S.TBDIS),'COVMATNorm',MultiCMNorm,...
        'fixPar','1 5 6',...
        'exclDataStart',exclDataStart);
    par(:,i) = F{i}.RESULTS{1};
    par(1,i) = par(1,i)+M.mnuSq_i;
    err(:,i) = F{i}.RESULTS{2};
    chi2min(i) = F{i}.RESULTS{3};
    dof(i) = F{i}.RESULTS{5};
    TBDIS_Fit(:,i) = F{i}.SO.TBDIS;
    if  abs(chi2min(i)-2.7)<ScanPrcsn % exit condition 
        fprintf('90 C.L. limit reached \n');
        mnuSq_i_Fit(i+1:end)=NaN; % set all not fitted to NaN
        chi2min(i+1:end) = NaN;
        dof(i+1:end)     = NaN;
        par(:,i+1:end)   = NaN;
        err(:,i+1:end)   = NaN;
        fprintf('Exist Scan, 90%% C.L. reached \n');
        break
    elseif chi2min(i)<2.7
        mNuSqLarger= min(mnuSq_i_Fit(mnuSq_i_Fit>mnuSq_i_Fit(i))); % find next higher mNuSq
        mnuSq_i_Fit(i+1) =  0.5*(mnuSq_i_Fit(i)+mNuSqLarger); % find middle between current and next higher one
    elseif chi2min(i)>2.7
        mNuSqSmaller= max(mnuSq_i_Fit(mnuSq_i_Fit<mnuSq_i_Fit(i))); % find next smaller mNuSq
        mnuSq_i_Fit(i+1) =  0.5*(mnuSq_i_Fit(i)+mNuSqSmaller);  % find middle between current and next smaller one
    end
end

% Compute Sensitivity
mNu90 =interp1(chi2min(~isnan(chi2min)),mnuSq_i_Fit(~isnan(chi2min)),min(chi2min(~isnan(chi2min)))+2.7,'spline'); % sensitivity on mnu^2 (90% C.L.)
mNumin = interp1(chi2min(~isnan(chi2min)),mnuSq_i_Fit(~isnan(chi2min)),min(chi2min(~isnan(chi2min))),'spline');   % mass with minimal chi2
toc
mNuPlot = 0:0.01:mNuStop^2;
chi2func = interp1(mnuSq_i_Fit(~isnan(chi2min)),chi2min(~isnan(chi2min)),mNuPlot,'spline');
%% plot
if strcmp(plotFit,'ON')
    close all;
    f11 = figure('Name','NuMassScanChi2Curve','Renderer','opengl');
    set(f11, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);    
    x = [0,mNu90];
    y = [min(chi2min),min(chi2min)+2.7]; 
    pInterp = plot(mNuPlot,chi2func,'-','LineWidth',3,'Color',rgb('CadetBlue'));
    hold on;
    pLine = plot(x,(min(chi2min)+2.7).*ones(numel(x),1),'--','LineWidth',3,'Color',rgb('IndianRed'));
    pScan = plot(mnuSq_i_Fit,chi2min,'o','MarkerSize',10,'MarkerFaceColor',rgb('Navy'),'MarkerEdgeColor',rgb('Navy'));
    plot(mNu90.*ones(numel(y),1),y,'--','LineWidth',3,'Color',rgb('IndianRed'));
    xlabel('m_{\nu}^2 (eV^2)')
    ylabel(['\chi^2 (',num2str(dof(1)),' dof)']);
    scan_leg = ['m_{\nu}^2  = ',sprintf('%.3f \\pm %.3f eV^2 (90%% C.L.)\n',mNumin,mNu90),...
        'm_{\nu}  = ',sprintf('%.3f \\pm %.3f eV   (90%% C.L.)',sqrt(mNumin),sqrt(mNu90))];
    line_leg ='\chi^2_{min}+2.7';
    leg = legend([pInterp, pLine],scan_leg,line_leg,'Location','northwest'); legend boxoff;
    if ~strcmp(chi2,'chi2Stat')
        leg.Title.String = ['stat + ',sprintf(strrep(SysEffect,'_','-'))];
    else
        leg.Title.String = 'stat';
    end
    PrettyFigureFormat; 
    title(sprintf('Samak Fit to Asimov \n - KATRIN %.0f years - MTD: %s - BKG %.0f mcps - B_{T},B_{max} = %.0f %%',...
        TimeSec/(60*60*24*365),strrep(TD,'_',' '),BKG_RateSec*1e3, WGTS_B_T/3.6*100));
    set(gca,'FontSize',18);
    xlim([0 mNu90*1.5]);
    ylim([0,interp1(mnuSq_i_Fit(~isnan(chi2min)),chi2min(~isnan(chi2min)),mNu90*1.5,'spline')]);
    save_name = sprintf('Sensitivity_NuMassScan_%s%s_%s_%.2fT-BT_%.2fyears_SysBudget%s',chi2,strrep(SysEffect,'_','-'),TD,WGTS_B_T,TimeSec/(60*60*24*365),SysBudget);
    if ~exist('../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/','dir')
        mkdir ../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/;
        mkdir ../sensitivity_nominalKATRIN/plots/pdf/NuMassScanChi2Curve/;
        mkdir ../sensitivity_nominalKATRIN/plots/fig/NuMassScanChi2Curve/;
    end
    export_fig(f11,['../sensitivity_nominalKATRIN/plots/png/NuMassScanChi2Curve/',save_name,'.png']);
    savefig(f11,['../sensitivity_nominalKATRIN/plots/fig/NuMassScanChi2Curve/',save_name,'.fig'],'compact');
    publish_figurePDF(f11,['../sensitivity_nominalKATRIN/plots/pdf/NuMassScanChi2Curve/',save_name,'.pdf']);
end
elseif strcmp(Scan,'OFF')
    F = FITC('SO',M,'DATA',Data,'fitter',fitter,...
        'chi2name',chi2,...
        'COVMAT', MultiCM+diag(S.TBDIS),'COVMATFrac', MultiCMFrac+diag(1./S.TBDIS),...
        'COVMATShape', MultiCMShape+diag(S.TBDIS),'COVMATNorm',MultiCMNorm,...
        'fixPar','5 6',...
        'exclDataStart',exclDataStart);
    par = squeeze(repmat(F.RESULTS{1},1,1,nFitMax));
    par(1,:) = par(1,:)+M.mnuSq_i;
    err = squeeze(repmat(F.RESULTS{2},1,1,nFitMax));
    chi2min = squeeze(repmat(F.RESULTS{3},1,1,nFitMax));
    dof = squeeze(repmat(F.RESULTS{5},1,1,nFitMax));
    mnuSq_i_Fit          = zeros(nFitMax,1);
    mNu90 = err(1)*1.64; % sensitivity on mnu^2 (90% C.L.) central, two-sided interval
    mNumin = par(1);
end
% OLD NuMassScan (slower)
%Init
% F = cell(length(mnuSq_i_Fit),1);
% par     = zeros(6,length(mnuSq_i_Fit));
% err     = zeros(6,length(mnuSq_i_Fit));
% chi2min = zeros(length(mnuSq_i_Fit),1);
% dof     = zeros(length(mnuSq_i_Fit),1);
% Data = [S.qU, S.TBDIS, sqrt(S.TBDIS)];
% % Do Scan
% if ischar(SysEffect)
% progressbar(sprintf('Nu-mass Scan %s %s',chi2,SysEffect));
% end
% for i =1:length(mnuSq_i_Fit)
%     progressbar(i/length(mnuSq_i_Fit));
%     M.mnuSq_i = mnuSq_i_Fit(i);
%     F{i} = FITC('SO',M,'DATA',Data,'fitter',fitter,...
%         'chi2name',chi2,...
%         'COVMAT', MultiCM+diag(S.TBDIS),'COVMATFrac', MultiCMFrac+diag(1./S.TBDIS),...
%         'COVMATShape', MultiCMShape,'COVMATNorm',MultiCMNorm,...
%         'fixPar','1 5 6',...
%         'exclDataStart',exclDataStart);
%     par(:,i) = F{i}.RESULTS{1};
%     par(1,i) = par(1,i)+M.mnuSq_i;
%     err(:,i) = F{i}.RESULTS{2};
%     chi2min(i) = F{i}.RESULTS{3};
%     dof(i) = F{i}.RESULTS{5};
%     TBDIS_Fit(:,i) = F{i}.SO.TBDIS; % Model Spectrum after the Fit
%     if strcmp(CutScanRange,'ON')
%         if chi2min(i) >= 3
%             mnuSq_i_Fit(i+1:end)=NaN; % set all not fitted to NaN
%             chi2min(i+1:end) = NaN; 
%             dof(i+1:end)     = NaN;
%             par(:,i+1:end)   = NaN;
%             err(:,i+1:end)   = NaN;
%             fprintf('Exist Scan, 90%% C.L. reached');
%             break
%         end
%     end
% end
%     par(2,:) = par(2,:)+M.Q_i;
%     par(4,:) = par(4,:)+M.BKG_RateSec_i;