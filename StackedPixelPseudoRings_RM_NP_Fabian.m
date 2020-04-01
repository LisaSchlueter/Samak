QAplots = 'OFF';
HVDrift = 'OFF';
Detrend = 'ON';


%% Rate Evolution --> mV equivalent
myMainTitle = sprintf('KATRIN - KNM2 - FPD @E_0-300eV - Non-Poissonian Components');
maintitle   = myMainTitle;
fig1      = figure('Name',sprintf('KATRIN - KNM2 Scanwise Background'),...
    'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

savefile1    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation.0f_PseudoRings_NP.png');


% Loop on Period
for j=1:3
    
    
    RunList   = ['KNM2_RW' num2str(j)];
    
    % Read Data
    DataType  = 'Real';
    RunAnaArg = {'RunList',RunList,'DataType',DataType,...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1};
    MR        = MultiRunAnalysis(RunAnaArg{:});
    range = 40;               % fit range in eV below endpoint        
    MR.exclDataStart = MR.GetexclDataStart(range); % find correct data, where to cut spectrum
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
    R         = A.MultiObj(1);
    R.ROIFlag='14keV'; R.SetROI;

    % Slow Control Data
    p1 = (R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2');
    p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);
    
    %% HV Drift Correction
    FirstDayPeriod1 = datetime('02-Oct-2019 14:52:19');
    TimeLineDaysFirstDayPeriod1 = days(MR.SingleRunData.StartTimeStamp-FirstDayPeriod1);
    if strcmp(HVDrift,'ON')
        HVdriftPerPixel = 1.5*(TimeLineDaysFirstDayPeriod1) * 6.3e-3;
    else
        HVdriftPerPixel = 0;
    end
    
    %% Stacked Pixel Data for each patch
    count = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateE = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    sstime    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    cf    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    for i=1:A.nRings
          R           = A.MultiObj(i);
          R.ROIFlag='14keV'; R.SetROI;
          R.RMCorrection('QAplots',QAplots);
          count  = R.SingleRunData.TBDIS_RM;
          sstime = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
          rate   = count./sstime;
          if strcmp(Detrend,'ON')
              rate = detrend(rate)+mean(rate);
          end
          corrcount_norm{j,i} = rate .* mean(sstime);
    
        if j==1
          subplot(3,4,i+j-1);
        elseif j==2
          subplot(3,4,4+i);
        elseif j==3
          subplot(3,4,8+i);
        end
        
        %% Stacked Pixel: Non-Poissonian Component?
        pdG = fitdist(corrcount_norm{j,i}','Normal');
        pdN = fitdist(corrcount_norm{j,i}','poisson');
        Broadening = sqrt(pdG.sigma^2-pdN.lambda)/(mean(sstime)*737.8) * 1e3 * 117 / numel(R.PixList);
        
        
        h = histogram(corrcount_norm{j,i},15,'Normalization','pdf',...
            'FaceColor',rgb('DodgerBlue'),'LineWidth',2,'FaceAlpha',0.7);
        xlabel(sprintf('counts in %.2f sec',mean(sstime)));
        ylabel('Frequency');
        PrettyFigureFormat
        % Gaussian PDF
        pdfG = @(b) 1/sqrt(2*pi*pdG.sigma.^2) .* exp(-(b-pdG.mu).^2/(2*pdG.sigma.^2));
        % Poisson PDF
        pdfP = @(b) 1/sqrt(2*pi*pdN.lambda) .* exp(-(b-pdN.lambda).^2/(2*pdN.lambda));
        hold on
        %b = (h.BinEdges(1:end-1)+h.BinWidth/2);
        b = linspace(min(corrcount_norm{j,i}),max(corrcount_norm{j,i}),100);
        g=plot(b,pdfG(b),'Color',rgb('IndianRed'),'LineWidth',4);
        p=plot(b,pdfP(b),'Color',rgb('GoldenRod'),'LineWidth',4);
        leg=legend([h g p],...
            sprintf('%.0f scans - PSR %.0f',numel(R.RunList),i),...
            sprintf('Gaussian'),...
            sprintf('Poisson \n NP = %.1f \n \\sigma_{eq} = %.2f mV',pdG.sigma/sqrt(pdN.lambda),Broadening),...
            'location','northwest');
        leg.Color = 'none'; legend boxoff;
        hold off
        PrettyFigureFormat
        set(gca,'FontSize',12);
        disp(pdG.sigma/sqrt(pdN.lambda));
        
    end
    export_fig(fig1,savefile1);

end