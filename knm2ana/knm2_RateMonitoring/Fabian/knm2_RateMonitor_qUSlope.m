% ------------------------------------------------------------
% KNM2 rate monitor point (-300eV) analysis
% calculate rate - qU dependence with MC simulations
% ------------------------------------------------------------
RecomputeFlag = 'OFF';
% You need to have Run GetSamakPath in your Working Directory Before
%function [Slope_RateqU,Err_Slope_RateqU] = knm2_RateMonitor_qUSlope
savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
MakeDir(savedir);

% Fake KNM2 Data - RW2 Confirguration
% All KNM2 RW2 inputs are in ref_FakeRun_KNM2_CD84_RateMonitor
% (fake run with 20 data points +-1eV around E0-300eV)
savename = [savedir,'knm2_RateMonitoring_qUSlope.mat'];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    InitFile =  @ref_FakeRun_KNM2_CD84_RateMonitor;
    CommonArg = {'FakeInitFile',InitFile,...
        'RunNr',1,...% has no meaning
        'DataType','Fake',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'chi2','chi2Stat',...
        'RingMerge','Full',...
        'minuitOpt','min;migrad',...
        'NonPoissonScaleFactor',1,...
        'AnaFlag','StackPixel',...
        'fixPar','mNu E0 Bkg Norm',...
        'RingList',1:12};
    % set up RunAnalysis object:
    A = RunAnalysis(CommonArg{:});
    % set ring segmentation
    %  R = RingAnalysis('RunAnaObj',A,'RingList',1:4);
    
    % get -300 eV rate / error
    qUIndex = find((A.RunData.qU-18574)<-95); % all Data points around rate monitor point
    qU      = A.RunData.qU(qUIndex);
    Rate    = A.RunData.TBDIS(qUIndex)./(A.RunData.qUfrac(qUIndex).*A.RunData.TimeSec);
    TBDIS = A.RunData.TBDIS(qUIndex);
    RateErr = sqrt(Rate);
    
    % linear fit of rate around E0-300V
    [par, err, chi2min,dof] = linFit(qU-18574,Rate,RateErr);
    
    save(savename,'par','err','chi2min','dof','qUIndex','qU','Rate','RateErr','CommonArg','TBDIS');
end

%% Sanity plot: Uniform - All Pixels
MeanRate = wmean(Rate,RateErr);

f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);
%     x = linspace(min(qU),max(qU),1e3);
%     y = x.*par(1)+par(2);
x = linspace(min(qU),max(qU),100)-18574;
y = par(1).*x+par(2);
p1 = plot(x,y.*1e-04,'LineWidth',2,'Color',rgb('Black'));
hold on;
pMC = errorbar(qU-18574,Rate.*1e-04,RateErr.*1e-04,'.','LineWidth',2,'MarkerSize',20,'CapSize',0);
hold on;
pNone = plot(NaN,NaN,'Color',rgb('White'));
pNone2 = plot(NaN,NaN,'Color',rgb('White'));

xlabel('Retarding energy - 18574 (eV)');
ylabel(sprintf('{\\itR}_{300}^{ sim.} (10^4 \\times cps)'));
PrettyFigureFormat('FontSize',18);
leg = legend([pMC,p1,pNone2],sprintf('Simulation (KNM-2-like)'),...
    sprintf('Slope: \\alpha = %.2f cps/mV ,',par(1)/1e3),...
    sprintf('with \\langle{\\itR}_{300}^{ sim.}\\rangle = %.2f \\times 10^4 cps',1e-04.*MeanRate));% sprintf('relative slope: a_{rel} = %.2f kV^{-1}',par(1)/MeanRate*1e3),...
leg.Location = 'southwest';

PrettyLegendFormat(leg,'alpha',0.7);
leg.FontSize = get(gca,'FontSize');
Slope_RateqU = par(1)/MeanRate*1e3;
Err_Slope_RateqU = err(1)/MeanRate*1e3;
xlim([-301 -298.8])
ylim([6.47 6.71]);
%
pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');

export_fig(pltname);
%print(pltname,'-dpng','-r400');% export_fig(pltname);
