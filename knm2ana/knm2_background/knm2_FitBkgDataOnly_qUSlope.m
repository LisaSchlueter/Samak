% knm1. fit only sub-scans above endpoint
% fit linear slope
% mnusq uncertainty as a function of slope constraints

% settings
DataType = 'Real';
chi2 = 'chi2Stat';
NP = 1.112;
savedir = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];
savename = [savedir,sprintf('knm2_FitBkgDataOnly_%s_%s_NP%.3f.mat',DataType,chi2,NP)];

if exist(savename,'file')
    load(savename);
else
    MakeDir(savedir);
     SigmaSq =  0.0124+0.0025;
     BKG_PtSlope = 3*1e-06;
    % Init Model Object and covariance matrix object
    A = MultiRunAnalysis('RunList','KNM2_Prompt',...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg',...
        'NonPoissonScaleFactor',NP,...
        'SysBudget',40,...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AngularTFFlag','ON',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag','FSD',...
        'FSD_Sigma',sqrt(SigmaSq),...
         'BKG_PtSlope',BKG_PtSlope);
    
    %% get scan-steps above endpoint
    qU      = A.RunData.qU(end-4:end)-18574;
    Time    = A.ModelObj.qUfrac(end-4:end).*A.ModelObj.TimeSec;
    BKG     = A.RunData.TBDIS(end-4:end)./Time;
    BKG_err = (sqrt(BKG.*Time).*NP)./Time;
    
    %    parInit       = [mean(BKG), 0];
    %    tmparg = sprintf(['set pri -10; min ; migrad '],'');
    %    Args   = {parInit, [qU,BKG,BKG_err], '-c', tmparg};
    %    [par2, err2,chi2min2, ~] = fminuit('chi2BKG',Args{:});
    
    [par, err,chi2min, ~] = linFit(qU,BKG-BKG(1),BKG_err);
    save(savename,'qU','BKG','BKG_err','Time','par','err','chi2min')
    
end

f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);
x = linspace(min(qU),max(qU),10);
pfit = plot(x,1e3.*(par(1).*x+par(2)+BKG(1)),'LineWidth',2,'Color',rgb('DarkOrange'));
hold on;
e1 = errorbar(qU,BKG.*1e3,BKG_err.*1e3,'.','MarkerSize',20,'LineWidth',2,'CapSize',0,'Color',rgb('DeepSkyBlue'));

% convert from cps/eV to mucps/keV
% =>  cps/eV  = (1e3 mcps)/(1e-3 keV) = 10^6 mcps/keV
% =>  cps/eV  = (1e3 mcps)/(eV)       = 10^3 mcps/eV
ConvFac = 1e6;

leg = legend([e1,pfit],'KNM2 data',...
    sprintf('Linear fit: {\\its}_{qU} = (%.0f\\pm%.0f) mcps/keV',par(1).*ConvFac,err(1).*ConvFac),...
    'Location','southeast');
PrettyLegendFormat(leg);

ylabel('Counts (mcps)');
xlabel('Retarding energy - 18574 (eV)');
xlim([-10 150])
ylim([218.2 225.5]);
PrettyFigureFormat('FontSize',20);
leg.FontSize = get(gca,'FontSize')+2;
%
pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = [pltdir,sprintf('knm2_FitBkgDataOnly_qUSlope.pdf')];
export_fig(pltname);
%%
fprintf('Slope significance %.1f sigma\n',abs(par(1)./err(1)));

