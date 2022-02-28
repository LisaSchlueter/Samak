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
    
    f1 = GetFigure;
    %     x = linspace(min(qU),max(qU),1e3);
    %     y = x.*par(1)+par(2);
    x = linspace(min(qU),max(qU),100)-18574;
    y = par(1).*x+par(2);
    p1 = plot(x,y.*1e-04,'LineWidth',2,'Color',rgb('Gray'));
    hold on;
    pMC = plot(qU-18574,Rate.*1e-04,'.','LineWidth',3,'MarkerSize',20);
    hold on;
    pNone = plot(NaN,NaN,'Color',rgb('White'));
    pNone2 = plot(NaN,NaN,'Color',rgb('White'));
    
    
    PrettyFigureFormat('FontSize',24);
    xlabel('Retarding energy - 18574 (eV)');
    ylabel(sprintf('Rate (10^4 cps)'));
%     leg = legend([pMC,p1],sprintf('MC simulation'),...
%         sprintf('Linear fit slope = %.2f \\pm %.2f cps/mV  = %.2f \\pm %.2f kV^{-1}',par(1)/1e3,err(1)/1e3,par(1)/MeanRate*1e3,err(1)/MeanRate*1e3));
leg = legend([pMC,p1,pNone,pNone2],sprintf('Simulation (KNM-2-like)'),...
        sprintf('Linear fit slope: a = %.2f cps\\cdotmV^{-1} ,',par(1)/1e3),...
        sprintf('relative slope: a_{rel} = %.2f kV^{-1}',par(1)/MeanRate*1e3),...
        sprintf('with \\langleRate @ 300eV\\rangle = %.2f\\cdot10^4 cps',1e-04.*MeanRate));
  leg.Location = 'southwest';

   PrettyLegendFormat(leg,'alpha',0.7);
    Slope_RateqU = par(1)/MeanRate*1e3;
    Err_Slope_RateqU = err(1)/MeanRate*1e3;
    xlim([-300.5 -299])
   
    pltdir = strrep(savedir,'results','plots');
    MakeDir(pltdir);
    pltname = strrep(strrep(savename,'results','plots'),'.mat','.png');
    
   print(pltname,'-dpng','-r400');% export_fig(pltname);
    