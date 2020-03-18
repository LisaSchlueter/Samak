function [PlasmaDrifts, Err_PlasmaDrifts, Offsets, Err_Offsets] = knm2_RingwiseFit(varargin)
    p = inputParser;
    p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('QAplots','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('ROI','14keV',@(x)ismember(x,{'14keV','Default'}));
    p.parse(varargin{:});
    saveplot = p.Results.saveplot;
    QAplots  = p.Results.QAplots;
    ROI      = p.Results.ROI;
    [zerolevel,Err_zerolevel] = ZeroLevel(ROI);
    [Slope_RateqU,Err_Slope_RateqU] = knm2_RateMonitor_qUSlope;
    Slope_RateqU = Slope_RateqU*1e-6;
    Err_Slope_RateqU = Err_Slope_RateqU*1e-6;
    OffsetMSCorr=[0, 29.1, 60, 91.6];
    for j = 1:3
        RunList   = ['KNM2_RW' num2str(j)];
        Knm2AnaArg = {'RunList',RunList,'DataType','Real',...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1,...
            'MosCorrFlag','OFF',...
            'TwinBias_Q',18573.7,...
            'ROIFlag',ROI};
        DataUni_RWn   = MultiRunAnalysis(Knm2AnaArg{:});
        range = 40;               % fit range in eV below endpoint        
        DataUni_RWn.exclDataStart = DataUni_RWn.GetexclDataStart(range); % find correct data, where to cut spectrum
        DataPSR_RWn   = RingAnalysis('RunAnaObj',DataUni_RWn,'RingList',1:4);

        for i = 1:DataUni_RWn.nRings
            DataPSR_RWn.MultiObj(i).RMCorrection('saveplot',saveplot,'pixlist',sprintf('ring%i',i),'QAplots',QAplots);
            [~,par,err,~,~] = DataPSR_RWn.MultiObj(i).PlotFitRunListCorr('Parameterx','time','Parametery','rate300','Fit','ON','saveplot',saveplot,'pixlist',sprintf('ring%i',i));
            MeanRate = mean(DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM./(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec))*numel(DataPSR_RWn.MultiObj(i).PixList)/numel(DataPSR_RWn.MultiObj(1).PixList);
            Err_MeanRate = std(DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM./(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec));
            Slope_RateTime = par(1)/MeanRate;
            Err_Slope_RateTime = err(1)/MeanRate;
            Offsets(i,j) = (MeanRate - zerolevel)/(Slope_RateqU*zerolevel)-OffsetMSCorr(i);
            RateOffsetError = sqrt(Err_MeanRate^2+Err_zerolevel^2);
            SlopeOffsetError = sqrt(Err_Slope_RateqU^2*zerolevel^2+Err_zerolevel^2*Slope_RateqU^2);
            Err_Offsets(i,j) = sqrt(RateOffsetError^2*(1/(Slope_RateqU*zerolevel))^2+SlopeOffsetError^2*((MeanRate - zerolevel)/((Slope_RateqU*zerolevel)^2))^2);
            PlasmaDrifts(i,j) = Slope_RateTime/Slope_RateqU*24;
            Err_PlasmaDrifts(i,j) = sqrt((Err_Slope_RateTime)^2*Slope_RateqU^-2+Err_Slope_RateqU^2*((Slope_RateTime)^2/(Slope_RateqU)^4))*24;
        end
    end
    fig88 = figure('Renderer','painters');
    set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
    PlotStyle = { '-o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),'LineWidth',2};
    e1 = errorbar(PlasmaDrifts,Err_PlasmaDrifts,PlotStyle{:});
    xlabel(sprintf('Pseudoring'));
    ylabel(sprintf('mV/day'));
    leg = legend('RW1','RW2','RW3');
    title(leg,sprintf('Plasma potential drifts'));
    leg.Location = 'best';
    PrettyFigureFormat;
    if strcmp(saveplot,'ON')
        SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/Knm2/PlasmaPotential/')];
        MakeDir(SaveDir);
        SaveName = sprintf('PlasmaDrifts_Ringwise.pdf');
        export_fig(fig88,[SaveDir,SaveName]);
    end
    fig88 = figure('Renderer','painters');
    set(fig88,'units','normalized','pos',[0.1, 0.1,1,0.6]);
    PlotStyle = { '-o','MarkerSize',8,'MarkerFaceColor',rgb('SkyBlue'),'LineWidth',2};
    e1 = errorbar(Offsets,Err_Offsets,PlotStyle{:});
    xlabel(sprintf('Pseudoring'));
    ylabel(sprintf('mV'));
    leg = legend('RW1','RW2','RW3');
    title(leg,sprintf('Offsets from RW1, Ring1'));
    leg.Location = 'best';
    PrettyFigureFormat;
    if strcmp(saveplot,'ON')
        SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/Knm2/PlasmaPotential/')];
        MakeDir(SaveDir);
        SaveName = sprintf('PlasmaPotential_Offsets_Ringwise.pdf');
        export_fig(fig88,[SaveDir,SaveName]);
    end
end

function [zerolevel,Err_zerolevel] = ZeroLevel(ROI)
    Knm2AnaArg = {'RunList','KNM2_RW1','DataType','Real',...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1,...
            'MosCorrFlag','OFF',...
            'TwinBias_Q',18573.7,...
            'ROIFlag',ROI};
        DataUni_RW1   = MultiRunAnalysis(Knm2AnaArg{:});
        range = 40;               % fit range in eV below endpoint        
        DataUni_RW1.exclDataStart = DataUni_RW1.GetexclDataStart(range); % find correct data, where to cut spectrum
        DataPSR_RW   = RingAnalysis('RunAnaObj',DataUni_RW1,'RingList',1);
        DataPSR_RW.MultiObj(1).RMCorrection('QAplots','OFF');
        zerolevel = mean(DataPSR_RW.MultiObj(1).SingleRunData.TBDIS_RM./(DataPSR_RW.MultiObj(1).SingleRunData.qUfrac_RM.*DataPSR_RW.MultiObj(1).SingleRunData.TimeSec));
        Err_zerolevel = std((DataPSR_RW.MultiObj(1).SingleRunData.TBDIS_RM./(DataPSR_RW.MultiObj(1).SingleRunData.qUfrac_RM.*DataPSR_RW.MultiObj(1).SingleRunData.TimeSec)));
end