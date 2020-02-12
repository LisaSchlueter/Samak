function [mnuSq_i_Fit, chi2minScan, dofScan,mNu90,mNumin] = OLD_NuMassScan_SensitivityNominal(varargin)
% ------------------------------------------------------------------------------------------------
% Goal: 
% Obtain sensitivity on neutrino mass squared
% Method: Scan
% Simulate KATRIN 1 time (asimov)
% Options: MTD, Background, Time
% Fit with statistic only or with a covariance matrix (E0,N,B free) and a fixed neutrino mass (scan over mnu)
% Plot chi^2 as a function of neutrino mass
% Compute sensitivity to neutrino mass squared
% ------------------------------------------------------------------------------------------------
tic
addpath(genpath('../../../Samak2.0'));
p = inputParser;
p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x) && x>0);
p.addParameter('TD','Sensitivity_30eV_Ba9G',@(x) ischar(x));
%Simulation Model Parameter
p.addParameter('MACE_Ba_T',9e-04,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_B_T',3.6*0.7,@(x)isfloat(x) && x>=0);
%Fit Model Paramater
p.addParameter('mnuSq_i_Fit',(0:20:200)*1e-03,@(x)all(isfloat(x))); %Scan range
p.addParameter('CutScanRange','ON',@(x)ismember(x,{'ON','OFF'}));  % Stop after chi2>=chi2+3
%Fit Settings
p.addParameter('chi2','chi2CM',@(x)ismember(x,{'chi2Stat','chi2CM'}));
p.addParameter('SysEffect','RF',@(x)ischar(x)); %1 sys effect e.g. FSD
p.addParameter('plotFit','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

TimeSec      = p.Results.TimeSec;
TD           = p.Results.TD;
MACE_Ba_T    = p.Results.MACE_Ba_T;
WGTS_B_T     = p.Results.WGTS_B_T;
mnuSq_i_Fit  = p.Results.mnuSq_i_Fit;
CutScanRange = p.Results.CutScanRange;
chi2         = p.Results.chi2;
SysEffect    = p.Results.SysEffect;
plotFit      = p.Results.plotFit;
if strcmp(chi2,'chi2Stat')
    SysEffect = '';
end
%---------------------------------------------parser end -----------------------------------------------%
%Init
chi2minScan = zeros(length(mnuSq_i_Fit),1);
dofScan =  zeros(length(mnuSq_i_Fit),1);

% Do Scan
progressbar(sprintf('Nu-mass Scan %s',SysEffect));
for i =1:length(mnuSq_i_Fit)
   progressbar(i/length(mnuSq_i_Fit));
    [par, err, chi2min, dof,qU,TBDIS_Sim,TBDIS_Fit, ModelObj] = ...
        SensitivityStudy_NominalKATRIN('StatFluct','OFF','SysFluct','OFF','nSamples',1,...
        'TD',TD,'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'chi2',chi2,'SysEffect',SysEffect,'TimeSec',TimeSec,...
        'mnuSq_i_Fit',mnuSq_i_Fit(i),'fixPar','1 5 6','saveResults','OFF','PlotFit','OFF');
    chi2minScan(i) = chi2min;
    dofScan(i) = dof;    
    if strcmp(CutScanRange,'ON')
        if chi2min >= 3
            mnuSq_i_Fit = mnuSq_i_Fit(1:i);
            chi2minScan = chi2minScan(1:i);
            dofScan = dofScan(1:i);
            fprintf('Exist Scan, 90%% C.L. reached');
            break
        end
    end
end

% Compute Sensitivity
mNu90 =interp1(chi2minScan,mnuSq_i_Fit,min(chi2minScan)+2.7,'spline'); % sensitivity on mnu^2 (90% C.L.)
mNumin = interp1(chi2minScan,mnuSq_i_Fit,min(chi2minScan),'spline');   % mass with minimal chi2
toc
%% plot
if strcmp(plotFit,'ON')
    close all;
    f11 = figure(11);
    set(f11, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
    
    x = [0,mNu90];
    y = [min(chi2minScan),min(chi2minScan)+2.7];
    pScan = plot(mnuSq_i_Fit,chi2minScan,'LineWidth',3,'Color',rgb('CadetBlue'));
    hold on;
    pLine = plot(x,(min(chi2minScan)+2.7).*ones(numel(x),1),'--','LineWidth',3,'Color',rgb('IndianRed'));
    plot(mNu90.*ones(numel(y),1),y,'--','LineWidth',3,'Color',rgb('IndianRed'));
    xlabel('m_{\nu}^2 (eV^2)')
    ylabel(['\chi^2 (',num2str(dofScan(1)),' dof)']);
    scan_leg = ['m_{\nu}^2 = ',sprintf('%.3f \\pm %.3f eV^2 (90%% C.L.)\n',mNumin,mNu90),...
        'm_{\nu}  = ',sprintf('%.3f \\pm %.3f eV (90%% C.L.)',sqrt(mNumin),sqrt(mNu90))];
    line_leg ='\chi^2_{min}+2.7';
    legend([pScan, pLine],scan_leg,line_leg,'Location','northwest'); legend boxoff;
    PrettyFigureFormat;
    set(gca,'FontSize',18);
    title(sprintf('Samak %s %s Fit to Simulation (asimov) \n - KATRIN %.0f years - MTD: %s -',chi2,SysEffect,TimeSec/(60*60*24*365),TD));
    
    save_name = sprintf('Sensitivity_NuMassScan_%s%s_%s_%.0f-BKG_%.2fyears',chi2,SysEffect,TD,MACE_Ba_T*1e3,TimeSec/(60*60*24*365));
    export_fig(f11,['./plots/png/',save_name,'.png']);
    savefig(f11,['./plots/fig/',save_name,'.fig'],'compact');
    publish_figurePDF(f11,['./plots/pdf/',save_name,'.pdf']);
end
end