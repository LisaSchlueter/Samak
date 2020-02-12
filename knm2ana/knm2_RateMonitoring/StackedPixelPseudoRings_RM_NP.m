% KNM2 - -300V Rate Monitor FPD
% Stakced-pixel Per Pseudo-Ring
% Evolution of Rate
% Rate Correction
% Conversion of rate in Potential Fluctuation
%
% Last Modified: 04/02/2020
% T. Lasserre
%

% Loop on Period
for j=1:3
    
    RunList   = ['KNM2_RW' num2str(j)];
    
    % Read Data
    DataType  = 'Real';
    RunAnaArg = {'RunList',RunList,'DataType',DataType,...
        'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel','RingMerge','Full'};
    MR        = MultiRunAnalysis(RunAnaArg{:});
    A         = RingAnalysis('RunAnaObj',MR,'RingList',1:4);
    R         = A.MultiObj(1);
    
    % Slow Control Data
    p1 = (R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')./mean((R.SingleRunData.WGTS_MolFrac_TT'+0.5*R.SingleRunData.WGTS_MolFrac_HT'+0.5*R.SingleRunData.WGTS_MolFrac_DT')).*R.SingleRunData.WGTS_CD_MolPerCm2'./mean(R.SingleRunData.WGTS_CD_MolPerCm2');
    p2 = mean(R.SingleRunData.qU_RM,1); p2=p2-mean(p2);

    %% Stacked Pixel Data for each patch
    count = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rate  = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    rateE = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    sstime    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    cf    = zeros(A.nRings,numel(A.RunAnaObj.RunList));
    for i=1:A.nRings
        R           = A.MultiObj(i);
        count(i,:)  = R.SingleRunData.TBDIS_RM;
        sstime(i,:) = mean(R.SingleRunData.qUfrac_RM,1).*R.SingleRunData.TimeSec;
        rate(i,:)   = count(i,:)./sstime(i,:);
        rateE(i,:)  = sqrt(count(i,:))./sstime(i,:);
        cf(i,:)     = R.RMRateErosCorrectionqUActivity;
        rateEquivalentmV_E{j,i}   =  (rateE(i,:) ./737.8 *1e3 * 117 / numel(R.PixList));
       
        count_norm{j,i}     = rate(i,:) .* mean(sstime(i,:));
        corrcount_norm{j,i} = count_norm{j,i}  .* cf(i,:) ;
    end
    
    
    %% Rate Evolution --> mV equivalent
    myMainTitle = sprintf('KATRIN - %s - FPD @E_0-300eV - Non-Poissonian Components',RunList);
    maintitle   = myMainTitle;
    fig1      = figure('Name',sprintf('KATRIN - %s Scanwise Background',RunList),...
        'NumberTitle','off','rend','painters','pos',[10 10 1400 1000]);
    a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
    a.FontSize=24;a.FontWeight='bold';
    
            savefile1    = sprintf('plots/KNM2_RM300_EffectivePotentialFluctuation_RW%.0f_PseudoRings_NP.png',j);

    for i=1:A.nRings
        
        %% Stacked Pixel: Non-Poissonian Component?
        pdG = fitdist(corrcount_norm{j,i}','Normal');
        pdN = fitdist(corrcount_norm{j,i}','poisson');
        
        subplot(A.nRings/2,A.nRings/2,i)
        
        h = histogram(corrcount_norm{j,i},15,'Normalization','pdf',...
            'FaceColor',rgb('DodgerBlue'),'LineWidth',2,'FaceAlpha',0.7);
        xlabel(sprintf('counts in %.2f sec',mean(sstime(i,:))));
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
            sprintf('%.0f subscans - Pseudo-Ring %.0f',numel(R.RunList),i),...
            sprintf('Gaussian \\sigma=%.2f counts',...
            (pdG.sigma)),sprintf('Poisson \\sigma=%.2f counts \n Non-Poisson Factor = %.1f',sqrt(pdN.lambda),pdG.sigma/sqrt(pdN.lambda)),...
            'location','northwest');
        leg.Color = 'none'; legend boxoff;
        hold off
        PrettyFigureFormat
        set(gca,'FontSize',12);
        disp(pdG.sigma/sqrt(pdN.lambda));
        
    end
    export_fig(fig1,savefile1);

end
