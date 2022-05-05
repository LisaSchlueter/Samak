% look at sensitivity contour from KNM-2 like simulation with flat MTD

for i=1:2
    if i==1
        FakeInitFile = @ref_KNM2_KATRIN_IsoStatMTD_new;
    else
        FakeInitFile = @ref_KNM2_KATRIN_IsoStatMTD_Bkg0mcps_new;
    end
    range = 40;
    freePar = 'E0 Norm Bkg';
    nGridSteps = 25;
    
    savedir = [getenv('SamakPath'),'ksn2ana/ksn2_MTD/results/'];
    FakeFileName = extractAfter(func2str(FakeInitFile),'ref_');
    savename = sprintf('%sContour_%s_%s.mat',savedir,strrep(freePar,' ',''),FakeFileName);
    
    if exist(savename,'file')
        load(savename)
    else
        % tritium run model
        F = RunAnalysis('RunNr',1,...
            'DataType','Fake',...
            'FakeInitFile',FakeInitFile,...
            'fixPar',freePar,...
            'SysBudget',40,...
            'AnaFlag','StackPixel',...
            'chi2','chi2Stat',...
            'FSDFlag','Sibille0p5eV',...
            'ELossFlag','KatrinT2A20',...
            'DopplerEffectFlag','FSD',...
            'RadiativeFlag','ON',...
            'minuitOpt','min ; minos');
        
        F.exclDataStart = F.GetexclDataStart(range);
        
        
        %% configure Sterile analysis object
        SterileArg = {'RunAnaObj',F,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
            'nGridSteps',nGridSteps,...
            'SmartGrid','OFF',...
            'RecomputeFlag','OFF',...
            'SysEffect','all',...
            'RandMC','OFF',...
            'range',range};
        S = SterileAnalysis(SterileArg{:});
        
        S.GridSearch('mNu4SqTestGrid',2);
        %     S.LoadGridFile('mNu4SqTestGrid',2);
        %     S.Interp1Grid;
        %    S.ContourPlot
        %     sin2T4_contour = S.sin2T4_contour;
        %     mNu4Sq_contour = S.mNu4Sq_contour;
        %     MakeDir(savedir);
        %     save(savename,'sin2T4_contour','mNu4Sq_contour');
    end
    
end

% GetFigure;
% plot(sin2T4_contour,mNu4Sq_contour);
% set(gca,'YScale','log');
% set(gca,'XScale','log');
