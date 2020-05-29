function [PlasmaDrifts, Err_PlasmaDrifts, Offsets, Err_Offsets] = knm2_RingwiseFit(varargin)
    p = inputParser;
    p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('QAplots','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('ROI','Default',@(x)ismember(x,{'14keV','Default'}));
    p.parse(varargin{:});
    saveplot = p.Results.saveplot;
    QAplots  = p.Results.QAplots;
    ROI      = p.Results.ROI;
    [ReferenceRate,~,ActivityRW2,qURW2] = ZeroLevel(ROI);
    %[Slope_RateqU,Err_Slope_RateqU] = knm2_RateMonitor_qUSlope %broken atm
    Slope_RateqU = -11.16*1e-6;        %cps/(eV*1 pixel)
    Err_Slope_RateqU = 1.43*1e-6;
    %Slope_RateqU = Slope_RateqU*1e-6;
    %Err_Slope_RateqU = Err_Slope_RateqU*1e-6;
    
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
        
        Activity{j} =  (DataPSR_RWn.MultiObj(1).SingleRunData.WGTS_MolFrac_TT+0.5*DataPSR_RWn.MultiObj(1).SingleRunData.WGTS_MolFrac_HT+0.5*DataPSR_RWn.MultiObj(1).SingleRunData.WGTS_MolFrac_DT)...
        .*(DataPSR_RWn.MultiObj(1).SingleRunData.WGTS_CD_MolPerCm2);

        for i = 1:DataUni_RWn.nRings
            
            DataPSR_RWn.MultiObj(i).RMCorrection('saveplot',saveplot,'pixlist',sprintf('ring%i',i),'QAplots',QAplots);
            Rate = DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM./(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec);
            RateCorrected = Rate./mean(Activity{j}).*ActivityRW2  + (mean(DataPSR_RWn.MultiObj(i).SingleRunData.qU_RM(1,:)) - qURW2(i)) * 6.3032 * numel(DataPSR_RWn.MultiObj(i).PixList);  % - ReferenceRate./numel(DataPSR_RWn.MultiObj(1).PixList).*numel(DataPSR_RWn.MultiObj(i).PixList);
            DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM = RateCorrected.*(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec);
            
            %[~,par,err,~,~] = DataPSR_RWn.MultiObj(i).PlotFitRunListCorr('Parameterx','time','Parametery','rate300','Fit','ON','saveplot',saveplot,'pixlist',sprintf('ring%i',i));
            Rate = DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM./(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec);
            Err_Rate = sqrt(DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM)./(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec);
            MeanRate = mean(DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM./(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec));
            Err_MeanRate = std(DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM./(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec));
            %Slope_RateTime = par(1)/MeanRate;
            %Err_Slope_RateTime = err(1)/MeanRate;
            DataPSR_RWn.MultiObj(i).SingleRunData.TBDIS_RM = -(Rate./737.8 * 1e3 * 117 / numel(DataPSR_RWn.MultiObj(i).PixList))...
                .*(DataPSR_RWn.MultiObj(i).SingleRunData.qUfrac_RM.*DataPSR_RWn.MultiObj(i).SingleRunData.TimeSec);
            [~,par,err,~,~] = DataPSR_RWn.MultiObj(i).PlotFitRunListCorr('Parameterx','time','Parametery','rate300','Fit','ON','saveplot',saveplot,'pixlist',sprintf('ring%i',i));
            NbxRunsPeriod = [121 95 92];
            
            Offsets(i,j) = mean(Rate./737.8 * 1e3 * 117 / numel(DataPSR_RWn.MultiObj(i).PixList));
            Err_Offsets(i,j) = mean(Err_Rate./737.8 * 1e3 * 117 / numel(DataPSR_RWn.MultiObj(i).PixList))./NbxRunsPeriod(j);
            PlasmaDrifts(i,j) = par(1)*24;
            Err_PlasmaDrifts(i,j) = err(1)*24;
            %PlasmaDrifts(i,j) = Slope_RateTime/Slope_RateqU*numel(DataPSR_RWn.MultiObj(i).PixList)/117*24;
            %Err_PlasmaDrifts(i,j) = sqrt((Err_Slope_RateTime)^2*Slope_RateqU^-2+Err_Slope_RateqU^2*((Slope_RateTime)^2/(Slope_RateqU)^4))*numel(DataPSR_RWn.MultiObj(i).PixList)/117*24;
        end
    end
    
    Offsets = -(Offsets-mean(mean(Offsets)));
    
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
    title(leg,sprintf('Offsets from baseline mean'));
    leg.Location = 'best';
    PrettyFigureFormat;
    if strcmp(saveplot,'ON')
        SaveDir = [getenv('SamakPath'),sprintf('tritium-data/plots/Knm2/PlasmaPotential/')];
        MakeDir(SaveDir);
        SaveName = sprintf('PlasmaPotential_Offsets_Ringwise.pdf');
        export_fig(fig88,[SaveDir,SaveName]);
    end
end

function [ReferenceRate,Err_ReferenceRate,ActivityRW2,qURW2] = ZeroLevel(ROI)
    Knm2AnaArg = {'RunList','KNM2_RW2','DataType','Real',...
            'FSDFlag','BlindingKNM2','ELossFlag','KatrinT2',...
            'AnaFlag','StackPixel','RingMerge','Full','NonPoissonScaleFactor',1,...
            'MosCorrFlag','OFF',...
            'TwinBias_Q',18573.7,...
            'ROIFlag',ROI};
        DataUni_RW2   = MultiRunAnalysis(Knm2AnaArg{:});
        range = 40;               % fit range in eV below endpoint        
        DataUni_RW2.exclDataStart = DataUni_RW2.GetexclDataStart(range); % find correct data, where to cut spectrum
        DataPSR_RW   = RingAnalysis('RunAnaObj',DataUni_RW2,'RingList',1);
        DataPSR_RW.MultiObj(1).RMCorrection('QAplots','OFF');
        qURW2 = zeros(DataUni_RW2.nRings,1);
        
        ReferenceRate = mean(DataPSR_RW.MultiObj(1).SingleRunData.TBDIS_RM./(DataPSR_RW.MultiObj(1).SingleRunData.qUfrac_RM.*DataPSR_RW.MultiObj(1).SingleRunData.TimeSec));
        Err_ReferenceRate = std((DataPSR_RW.MultiObj(1).SingleRunData.TBDIS_RM)./(DataPSR_RW.MultiObj(1).SingleRunData.qUfrac_RM.*DataPSR_RW.MultiObj(1).SingleRunData.TimeSec));
        ActivityRW2  =  mean(DataPSR_RW.MultiObj(1).SingleRunData.WGTS_MolFrac_TT'+0.5*DataPSR_RW.MultiObj(1).SingleRunData.WGTS_MolFrac_HT'...
            +0.5*DataPSR_RW.MultiObj(1).SingleRunData.WGTS_MolFrac_DT').*mean(DataPSR_RW.MultiObj(1).SingleRunData.WGTS_CD_MolPerCm2);
        
        for i = 1:DataUni_RW2.nRings
            qURW2(i) = mean(DataPSR_RW.MultiObj(i).SingleRunData.qU_RM(1,:));
        end
        
end