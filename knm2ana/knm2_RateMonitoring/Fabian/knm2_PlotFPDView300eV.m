% plot rate monitor point (KNM1: -200eV , KNM2: -300eV)
% FPD viewer
% Lisa, March 2020
RunList = 'Thierry';%'KNM1rm';
Corrections = 'OFF'; % Type of rate corrections
savedir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/results/'];
savename = sprintf('RateMonitorPoint_FPDView_%s_Corr%s.mat',RunList,Corrections);
Mode = 'Ring'; % 'Pixel'
if exist([savedir,savename],'file') && ~strcmp(RunList,'Thierry')
    load([savedir,savename])
else
    
    if strcmp(RunList,'Thierry')

        savenameT = sprintf('%sPixelWise_RM_%s_KNM1CorFlag%s_HVdriftCorFlag%s.mat',...
        savedir,'KNM1rm','OFF','OFF');
        d = importdata(savenameT);
        qU_RM = d.qU_RM;
        TBDIS_RM = NaN.*zeros(148,1);
        TBDIS_RM(d.PixList) = sum(PixelMap,2);
        RingPixList = d.RingPixList;
    else
        RunAnaArg = {'RunList',RunList,... % all KNM2 golden runs
            'DataType','Real',...
            'ROIFlag','14keV',...
            'NonPoissonScaleFactor',1};
        
        %% build object of MultiRunAnalysis class
        A = MultiRunAnalysis(RunAnaArg{:});
        A.ReadSingleRunData;
        TBDIS_RM =  zeros(148,1).*NaN;
        CDcorr = A.RunData.WGTS_CD_MolPerCm2./A.SingleRunData.WGTS_CD_MolPerCm2;
        switch Corrections
            case 'OFF'
                CorrFactor =  ones(numel(A.PixList),1);
            case 'CD'
                CorrFactor =  CDcorr;
        end
        
        Time_RM = A.RunData.qUfrac_RM.*A.RunData.TimeSec;
        TBDIS_RM(A.PixList) = sum(CorrFactor.*A.SingleRunData.TBDIS_RM(A.PixList,:),2)./Time_RM; % stack runs
        PixList = A.PixList;
        qU_RM = A.RunData.qU_RM;
        RingPixList = A.RingPixList;
        save(savename,'TBDIS_RM','Time_RM','CorrFactor','PixList','qU_RM','RingPixList')
    end
end

%%
meanTBDIS_RM = zeros(size(RingPixList,1),1);
switch Mode
    case 'Ring'
        for i=1:size(RingPixList,1)
            TBDIS_RM(RingPixList{i}) = mean(TBDIS_RM(RingPixList{i}));
            meanTBDIS_RM(i) = mean(TBDIS_RM(RingPixList{i}));
        end
    case 'Pixel'
end

meanTBDIS_RM = meanTBDIS_RM./nanmean(TBDIS_RM);
%% plot
fprintf('Rate at qU -18574 = %.0f eV \n',qU_RM-18574)
[plotHandle, cbHandle] = FPDViewer(TBDIS_RM./nanmean(TBDIS_RM));
cbHandle.Label.String = sprintf('Relative rate');
cbHandle.Label.FontSize = 24;
cbHandle.Ticks = 0.998:0.001:1.002;
cbHandle.Limits = ([0.9975,1.0022]);
%% save
plotdir = [getenv('SamakPath'),'knm2ana/knm2_RateMonitoring/plots/'];
MakeDir(plotdir);
plotname = [plotdir,strrep(savename,'.mat',[Mode,'.pdf'])];
fprintf('save plot to %s \n',plotname);
export_fig(plotHandle,plotname);
print(gcf,strrep(plotname,'.pdf','.png'),'-r400','-dpng');