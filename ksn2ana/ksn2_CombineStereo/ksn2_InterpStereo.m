% interpolate Stereo Chi^2 Map to match KATRIN grid
% sanity plots
InterpMode = 'spline';
Maxm4Sq = 40^2;
DataType = 'Real';
 SanityPlt = 'ON';

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefile = sprintf('%sksn2_%s_InterpStereo_%s_Max%.0feV2.mat',...
    savedir,DataType,InterpMode,Maxm4Sq);

if exist(savefile,'file') && 1==2
    fprintf('file exist %s \n',savefile);
else
   
    chi2 = 'chi2CMShape';
    
    %% settings that might change
    range = 40;
    if strcmp(DataType,'Real')
        nGridSteps = 40;%50;
        LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
    else
        nGridSteps = 40;
        LoadGridArg =   {'mNu4SqTestGrid',5,'ExtmNu4Sq','OFF'};
    end
    %% configure RunAnalysis object
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.112;
    end
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    %% configure Sterile analysis object
    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',LoadGridArg};
    %% load KATRIN (interpolated) grid
   
    S = SterileAnalysis(SterileArg{:});
    S.LoadGridFile(S.LoadGridArg{:});
    S.InterpMode = InterpMode;
    S.Interp1Grid('Maxm4Sq',Maxm4Sq);
    S.ContourPlot('BestFit','ON'); close;
    mNu4Sq_Katrin = S.mNu4Sq;
    sin2T4_Katrin = S.sin2T4;
    chi2_Katrin   = S.chi2-S.chi2_ref;
    chi2_ref_Katrin = S.chi2_ref;
    mNu4Sq_Katrin_bf = S.mNu4Sq_bf;  mNu4Sq_Katrin_contour = S.mNu4Sq_contour;
    sin2T4_Katrin_bf = S.sin2T4_bf;  sin2T4_Katrin_contour = S.sin2T4_contour;
    
    [~,sin2T4_Katrin_Osc] = S.Convert2Osci; % osci par space
    
    %% load Stereo
    StereoFile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/StereoMaps.mat'];
    ds = importdata(StereoFile);
    
    % find indices (in KATRIN grid) that match parameter space covered by STEREO
    IntermNuIdx = mNu4Sq_Katrin>=min(min(ds.mNu41Sq)) &  mNu4Sq_Katrin<=max(max(ds.mNu41Sq));
    InterSinIdx = sin2T4_Katrin_Osc>=min(min(ds.sin2TSq)) &  sin2T4_Katrin_Osc<=max(max(ds.sin2TSq));
    InterIdx = IntermNuIdx.* InterSinIdx;
    
    if strcmp(SanityPlt,'ON')
        % sanity plot: show which part of KATRIN parameter space is also
        % covered with STEREO
        GetFigure;
        surf(sin2T4_Katrin_Osc,mNu4Sq_Katrin,InterIdx,'EdgeColor','none')
        grid off
        view(2)
        set(gca,'Xscale','log')
        set(gca,'yscale','log')
        xlabel(sprintf('sin^2(2\\theta_{4})'));
        ylabel(sprintf('\\Deltam_{41}^2 (eV^2)'));
        PrettyFigureFormat;
        title('yellow = Stereo parameter space')
        ylim([0 40^2])
        xlim([0 1]);
        pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombinedStereo/plots/']; 
        MakeDir(pltdir);
        print(gcf,[pltdir,'StereoParSpace.png'],'-dpng','-r300');
        
    end
    
    %% interpolate STEREO
    IdxmNu = find(InterIdx(1000,:));%,1,'last');%find(InterIdx(:,1),1);
    Idxsin = find(InterIdx(:,10));%,1);
    mNu4Sq_tmp = mNu4Sq_Katrin(1,IdxmNu);%
    sin2T4_tmp = sin2T4_Katrin_Osc(Idxsin,1);% 
    mNu4Sq_inter = repmat(mNu4Sq_tmp,numel(sin2T4_tmp),1); % part of KATRIN that is interpolated
    sin2T4_inter = repmat(sin2T4_tmp,1,numel(mNu4Sq_tmp));  
    [X,Y] = meshgrid(ds.mNu41Sq(1,:),ds.sin2TSq(:,1));
    
    % stereo only
    chi2_inter = interp2(ds.mNu41Sq,ds.sin2TSq,ds.chi2,mNu4Sq_inter,sin2T4_inter,'lin');
    chi2crit_inter = interp2(ds.mNu41Sq,ds.sin2TSq,ds.chi2crit,mNu4Sq_inter,sin2T4_inter,'lin');

    %% show interpolation of STEREO
    
    if strcmp(SanityPlt,'ON')
        close
        figure('Units','normalized','Position',[0.1,0.1,0.8,0.5]);
        s1=  subplot(1,2,1);
        surf(ds.sin2TSq,ds.mNu41Sq,ds.chi2,'EdgeColor','none')
        grid off
        view(2)
        set(gca,'Xscale','log')
        set(gca,'yscale','log')
        xlabel(sprintf('sin^2(2\\theta_{4})'));
        ylabel(sprintf('\\Deltam_{41}^2 (eV^2)'));
        PrettyFigureFormat;
        c2 = colorbar;
        c2.Label.String = sprintf('\\Delta\\chi^2');
        c2.Label.FontSize = get(gca,'FontSize');
        title('STEREO original :  Phys.Rev.D 102 (2020) 052002','FontSize',15,'FontWeight','normal');
        
        s2 = subplot(1,2,2);
        surf(sin2T4_inter,mNu4Sq_inter,chi2_inter,'EdgeColor','none')
        grid off
        view(2)
        set(gca,'Xscale','log')
        set(gca,'yscale','log')
        xlabel(sprintf('sin^2(2\\theta_{4})'));
        ylabel(sprintf('\\Deltam_{41}^2 (eV^2)'));
        PrettyFigureFormat;
        c2 = colorbar;
        c2.Label.String = sprintf('\\Delta\\chi^2');
        c2.Label.FontSize = get(gca,'FontSize');
        title('STEREO interpolated','FontSize',15,'FontWeight','normal');
          
        linkaxes([s1,s2],'xy');
        pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombinedStereo/plots/'];
        MakeDir(pltdir);
        print(gcf,[pltdir,'StereoInterp.png'],'-dpng','-r300');
        
    end
    
    %% convert stero map to have same size as KATRIN
    % par space that is not covered in STEREO (but in KARTIN) has values=0
    chi2Stereo = zeros(size(S.chi2));
    chi2Stereo(logical(InterIdx)) = chi2_inter;
    % absolute values
    chi2StereoAbs = zeros(size(S.chi2));
    chi2StereoAbs(logical(InterIdx)) = chi2_inter+128.4 ; %(112dof)
    
    chi2critStereo = 5.99.*ones(size(S.chi2));
    chi2critStereo(logical(InterIdx)) = chi2crit_inter;
    
    if strcmp(SanityPlt,'ON')
        figure('Units','normalized','Position',[0.1,0.1,0.8,0.5]); 
        s1=  subplot(1,2,1);
        surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2_Katrin,'EdgeColor','none')
        hold on;
       pcontour =plot3(sin2T4_Katrin_contour,mNu4Sq_Katrin_contour,ones(numel(mNu4Sq_Katrin_contour),1).*1e4,...
         'LineWidth',2,'Color',rgb('White'));
       pbf = plot3(sin2T4_Katrin_bf,mNu4Sq_Katrin_bf,1e4,'wx','LineWidth',2);
        grid off
        view(2)
        set(gca,'Xscale','log')
        set(gca,'yscale','log')
        xlabel(sprintf('|{\\itU}_{e4}|^2'));
        ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
        PrettyFigureFormat;
        c1 = colorbar;
        c1.Label.String = sprintf('\\Delta\\chi^2');
        c1.Label.FontSize = get(gca,'FontSize');
        title('KATRIN','FontSize',15,'FontWeight','normal');
        leg = legend(pcontour,'KATRIN exclusion at 95% C.L.','Location','southwest');
        PrettyLegendFormat(leg);
        
        s2 = subplot(1,2,2);
        chi2StereoPlt = NaN.*zeros(size(S.chi2));
        chi2StereoPlt(logical(InterIdx)) = chi2_inter;
        surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2StereoPlt,'EdgeColor','none') %
        grid off
        view(2)
        set(gca,'Xscale','log')
        set(gca,'yscale','log')
        xlabel(sprintf('|{\\itU}_{e4}|^2'));
        ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
        PrettyFigureFormat;
        title('STEREO','FontSize',15,'FontWeight','normal');
        c2 = colorbar;
        c2.Label.String = sprintf('\\Delta\\chi^2');
        c2.Label.FontSize = get(gca,'FontSize');
        
        linkaxes([s1,s2],'xy');
        ylim([0 40^2]);
        xlim([1e-03 0.5])
        
        pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombinedStereo/plots/'];
        MakeDir(pltdir);
        print(gcf,[pltdir,'KATRINandStereo.png'],'-dpng','-r300');  
    end
    
    % combine
    chi2Combi    = chi2Stereo+chi2_Katrin;
    chi2CombiAbs = chi2StereoAbs+chi2_Katrin+S.chi2_ref;
    
    %% find minimum in shared parameter space
    InterIdx = logical(InterIdx);
    Chi2_refCombi =  min(min(chi2Combi(InterIdx)));
    chi2Combi(InterIdx) = chi2Combi(InterIdx)- Chi2_refCombi;
   
    %%
    save(savefile,'chi2Combi','chi2CombiAbs','Chi2_refCombi',...%combi
        'chi2Stereo','chi2critStereo','chi2StereoAbs',...%STEREO only
        'chi2_Katrin','chi2_ref_Katrin','mNu4Sq_Katrin','sin2T4_Katrin',...%KARTIN Only
        'mNu4Sq_Katrin_bf','mNu4Sq_Katrin_contour','sin2T4_Katrin_bf',...%KARTIN Only
        'sin2T4_Katrin_contour','sin2T4_Katrin_Osc',...%KARTIN Only
        'InterIdx');%information which part is covered in STEREO (in 1000x1000 matrix)
        
end


