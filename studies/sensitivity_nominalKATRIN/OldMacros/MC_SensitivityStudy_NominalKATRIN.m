function [par, err, chi2min, dof,qU,TBDIS_Sim,TBDIS_Fit, M] = MC_SensitivityStudy_NominalKATRIN(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity Test for KATRIN nominal settings
% Compute simulated spectra (amount: nSamples)
% - Options: with Stat and/or Sys Fluctuations
% Fit these spectra
% - Options: -with or without systematics (Covariance Matrix, different Effects or all)
%            -some Parameters fixed (1=neutrino mass, 2=endpoint,3=Background, 4 = Normalization, 5&6= FSD Parameters) 
% Output = Fit Results, Simulated Spectra, Fit Spectra, ModelObject from Fit            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=inputParser;
p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x) && x>0);
p.addParameter('TD','Flat60',@(x)ischar(x));
%Simulation Model Parameter
p.addParameter('MACE_Ba_T',3e-04,@(x)isfloat(x) && x>=0);
p.addParameter('WGTS_B_T',3.6*0.7,@(x)isfloat(x) && x>=0);
p.addParameter('StatFluct','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SysFluct','ON',@(x)ismember(x,{'ON','OFF'}));
%Fit Model Paramater
p.addParameter('mnuSq_i_Fit',0,@(x)isfloat(x));

%Fit Settings
p.addParameter('fixPar','5 6',@(x)ischar(x));
p.addParameter('chi2','chi2Stat',@(x)ismember(x,{'chi2Stat','chi2CM'}));
p.addParameter('SysEffect','FSD',@(x)ischar(x)); %1 sys effect e.g. FSD
p.addParameter('nSamples',100,@(x)isfloat(x));
p.addParameter('saveResults','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('plotFit','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

TimeSec     = p.Results.TimeSec;
TD          = p.Results.TD;
MACE_Ba_T   = p.Results.MACE_Ba_T;
WGTS_B_T    = p.Results.WGTS_B_T;
StatFluct   = p.Results.StatFluct;
SysFluct    = p.Results.SysFluct;
mnuSq_i_Fit = p.Results.mnuSq_i_Fit;
fixPar      = p.Results.fixPar;
chi2        = p.Results.chi2;
nSamples    = p.Results.nSamples;
saveResults = p.Results.saveResults;
plotFit     = p.Results.plotFit;
SysEffect   = p.Results.SysEffect;

%Background
BKG_RateSim = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T);

% Fit settings
fitter = 'minuit';
exclDataStart = 1;

if strcmp(chi2,'chi2Stat')
    SysEffect = '';  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%Simulation: TBD Object with nominal KATRIN settings
S = ref_TBD_NominalKATRIN('BKG_RateAllFPDSec',BKG_RateSim,'MACE_Ba_T',MACE_Ba_T,...
                          'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',(WGTS_B_T/3.6)*6,...
                           'TimeSec',TimeSec,'TD',TD);
S.ComputeTBDDS;
S.ComputeTBDIS;

%Model for Fit
M = ref_TBD_NominalKATRIN('BKG_RateAllFPDSec',BKG_RateSim,'MACE_Ba_T',MACE_Ba_T,...
                          'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',(WGTS_B_T/3.6)*6,...
                          'TimeSec',TimeSec,'TD',TD,...
                          'mnuSq_i',mnuSq_i_Fit);
M.ComputeTBDDS; M.ComputeTBDIS;
TBDIS_Fit = zeros(M.nqU,nSamples); %Save Fit spectra

%Covariance Matrix
if ~strcmp(chi2,'chi2Stat')
    [~, MultiCM, MultiCMFrac, MultiCMShape, MultiCMNorm] = ...
        ref_CovarianceMatrix_NominalKATRIN('ModelObj',M,'SysEffect',SysEffect);
elseif strcmp(chi2,'chi2Stat')
    MultiCM      = 0;
    MultiCMFrac  = 0;
    MultiCMShape = 0;
    MultiCMNorm  = 0;
end

%Compute MC Data set
switch StatFluct
    case 'ON'
        if strcmp(SysFluct,'OFF')
            TBDIS_Sim = mvnrnd(S.TBDIS,S.TBDIS',nSamples)'; % nSamples simulated integral spectra
        elseif strcmp(SysFluct,'ON')
            TBDIS_Sim = mvnrnd(S.TBDIS,MultiCM+diag(S.TBDIS),nSamples)';
        end
    case 'OFF'
        if strcmp(SysFluct,'OFF')
            TBDIS_Sim = repmat(S.TBDIS,1,nSamples);
        elseif strcmp(SysFluct,'ON')
            TBDIS_Sim = mvnrnd(S.TBDIS,MultiCM,nSamples)';
        end       
end

%Init variables
par     = zeros(6,nSamples);
err     = zeros(6,nSamples);
chi2min = zeros(nSamples,1);
dof     = zeros(nSamples,1);
F = cell(nSamples,1);

% Fit
qU = S.qU;
parfor i = 1:nSamples
Data = [qU, TBDIS_Sim(:,i), sqrt(TBDIS_Sim(:,i))];
F{i} = FITC('SO',M,'DATA',Data,'fitter',fitter,...
    'chi2name',chi2,...
    'COVMAT', MultiCM+diag(TBDIS_Sim(:,i)),'COVMATFrac', MultiCMFrac+diag(1./TBDIS_Sim(:,i)),...
    'COVMATShape', MultiCMShape,'COVMATNorm',MultiCMNorm,...
    'fixPar',fixPar,...
    'exclDataStart',exclDataStart);
par(:,i) = F{i}.RESULTS{1};
err(:,i) = F{i}.RESULTS{2};
chi2min(i) = F{i}.RESULTS{3};
dof(i) = F{i}.RESULTS{5};
TBDIS_Fit(:,i) = F{i}.SO.TBDIS; % Model Spectrum after the Fit
end

% save and plot

if strcmp(saveResults,'ON')
    save_name = sprintf('./results/SensitivityStudy_NominalKATRIN_MTD-%s_Time-%.2fy_mnu%.1feV_StatFluct%s_SysFluct%s_%.0fSamples_%s%s.mat',...
        TD,TimeSec/(60*60*24*365),mnuSq_i_Fit,StatFluct,SysFluct,nSamples,chi2,SysEffect);
    save(save_name,'par','err','chi2min','dof','qU','TBDIS_Sim','TBDIS_Fit');
end
if strcmp(plotFit,'ON')
    if nSamples>1
    else
    close all;
    f7 = figure(7);
    set(f7, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 1]);
    % TBDIS_Sim = TBDIS_Sim./(S.TimeSec*S.qUfrac);
    % TBDIS_Fit = TBDIS_Fit./(M.TimeSec.*M.qUfrac);
    s1 = subplot(2,1,1);
    pfit = plot(qU-18575,TBDIS_Fit,'-','LineWidth',3,'Color',rgb('CadetBlue'));
    hold on;
    pdata = errorbar(qU-18575,TBDIS_Sim,sqrt(TBDIS_Sim),'ko','MarkerSize',6,'MarkerFaceColor',.9*[1 1 1]);
    
    myleg = cell(length(par),1);
    i_factor = [1, 1, 1e3, 1];
    i_add    = [M.mnuSq_i, M.Q_i,M.BKG_RateSec_i,1];
    i_unit   = {'eV^2';'eV'; 'mcps';''};
    for i=1:length(par)
        if err(i)==0
            myleg{i}=sprintf('fixed');
        else
            myleg{i} =sprintf('= %.2f \t\\pm %.2f %s',(par(i)+i_add(i))*i_factor(i),err(i),i_unit{i});
        end
    end
    fit_leg = sprintf('Fit (%s %s):\n \\chi2 / dof=%.2f/%.0f \n m^2 %s \n E0_{eff} %s \n B %s \n N %s',...
        chi2,SysEffect,chi2min,dof,myleg{1},myleg{2},myleg{3},myleg{4});
    legend([pdata,pfit],'Simulated Data',fit_leg);
    xlabel('retarding potential - 18575 (V)');
    ylabel('counts')
    PrettyFigureFormat;
    set(gca,'FontSize',16);
    xlim([min(qU)-18575 max(qU)-18575]);
    title(sprintf('Samak Fit Simulation KATRIN %.0f years \n-MTD: %s - BKG %.0fmcps - %s',TimeSec/(60*60*24*365),TD,BKG_RateSim*1e3,StatFluctLabel));
    
   s2 = subplot(2,1,2);
   plot(qU-18575,(TBDIS_Sim-TBDIS_Fit)./sqrt(TBDIS_Sim),'ko','MarkerSize',6,'MarkerFaceColor',.9*[1 1 1]);
   hold on;
   plot(qU-18575,zeros(numel(qU),1),'k--');
   xlabel('retarding potential - 18575 (V)');
   ylabel('norm. residuals (\sigma)');
   PrettyFigureFormat;
   set(gca,'FontSize',16);
   xlim([min(qU)-18575 max(qU)-18575]);
   linkaxes([s1, s2],'x');
   
   save_name = sprintf('FitSimulation_%s_%.0f-BKG_%.2fyears_%s_%s%s',TD,BKG_RateSim*1e3,TimeSec/(60*60*24*365),StatFluctLabel,chi2,SysEffect);
   export_fig(f7,['./plots/png/',save_name,'.png']);
   savefig(f7,['./plots/fig/',save_name,'.fig'],'compact');
   publish_figurePDF(f7,['./plots/pdf/',save_name,'.pdf']);
    end
end
end