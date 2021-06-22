DataType = 'Real';%'Twin';
Maxm4Sq    = 38^2; % interpolation
freePar  = 'E0 Norm Bkg';
Plt = 'Contour';

savedir  = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/'];
savefile = sprintf('%sksn2_CombineSTEREO_%s_%s.mat',...
    savedir,DataType,strrep(freePar,' ',''));

if exist(savefile,'file')  && 1==2
    load(savefile);
else
    %% configure RunAnalysis object
    range = 40;
    chi2 = 'chi2CMShape';
    if strcmp(chi2,'chi2Stat')
        NonPoissonScaleFactor = 1;
    elseif  strcmp(chi2,'chi2CMShape')
        NonPoissonScaleFactor = 1.112;
    end
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...%free par
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
    %%
    nGridSteps = 40;
    switch DataType
        case 'Real'
            LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
            legStr = 'Data';
            BF = 'ON';
        case 'Twin'
            if ~contains(freePar,'mNu')
                nGridSteps = 30;
                LoadGridArg = {'mNu4SqTestGrid',2,'ExtmNu4Sq','ON'};
            else
                nGridSteps = 40;
                LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
            end
            legStr = 'MC Twin';
            BF = 'OFF';
    end
    
    %% configure Sterile analysis object
    SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
        'nGridSteps',nGridSteps,...
        'SmartGrid','OFF',...
        'RecomputeFlag','OFF',...
        'SysEffect','all',...
        'RandMC','OFF',...
        'range',range,...
        'LoadGridArg',LoadGridArg};
    S = SterileAnalysis(SterileArg{:});
    %% load KATRIN contour
    S.LoadGridFile(S.LoadGridArg{:});
    S.Interp1Grid('Maxm4Sq',Maxm4Sq);
    S.ContourPlot('BestFit',BF)
    close;
    sin2T4_contour_KATRIN = S.sin2T4_contour;
    mNu4Sq_contour_KATRIN = S.mNu4Sq_contour;
    if strcmp(BF,'ON')
        chi2_bf_KATRIN        = S.chi2_bf;
        sin2T4_bf_KATRIN      = S.sin2T4_bf;
        mNu4Sq_bf_KATRIN      = S.mNu4Sq_bf;
    end
    sin2T4 = S.sin2T4;
    mNu4Sq = S.mNu4Sq;
    chi2_KATRIN = S.chi2;
    %% load STEREO Minifile
    savedir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/'];
    if strcmp(DataType,'Real')
    savefileSTEREO = sprintf('%sksn2_InterpStereoMini_Max%.0feV2_Min%.2gfeV2.mat',...
        savedir,Maxm4Sq,min(min(mNu4Sq)));
    else
         savefileSTEREO = sprintf('%sksn2_InterpStereoMini_Max%.0feV2_Min%.2gfeV2_Sensitivity.mat',...
        savedir,Maxm4Sq,min(min(mNu4Sq)));
    end
    
    if exist(savefileSTEREO,'file')
        fprintf('load file %s \n',savefileSTEREO);
        dS = importdata(savefileSTEREO);
    else
        return
    end
    
    % double check : same binning?
    if any(any(dS.mNu4Sq~=mNu4Sq)) || any(any(dS.sin2T4~=sin2T4))
        fprintf('Stereo and KATRIN dont have same binning! \n');
        return
    end
    %% get STEREO only contour
    sin2T4_contour_STEREO = zeros(1e3,1);
    mNu4Sq_contour_STEREO = mNu4Sq(1,:);
    for i=1:1e3
        Idx = find(dS.chi2Stereo(:,i)>=dS.chi2Stereocrit(:,i),1,'first');
        if isempty(Idx)
            sin2T4_contour_STEREO(i) = NaN;
        else
            sin2T4_contour_STEREO(i) = sin2T4(Idx,i);
        end
    end
  
    if strcmp(DataType,'Real')
        chi2_bf_STEREO        = min(min(dS.chi2Stereo_cut));
        [row,col]             = find(dS.chi2Stereo_cut==chi2_bf_STEREO);
        sin2T4_bf_STEREO      = dS.sin2T4_cut(row,col);
        sin2T4_Osci_bf_STEREO = dS.sin2T4_Osci_cut(row,col);
        mNu4Sq_bf_STEREO      = dS.mNu4Sq_cut(row,col);
    end
    %% build shared chi^2 map
    if strcmp(DataType,'Real')
%         % 2 regions:
%         %1) KATRIN + STEREO
%         Chi2Stereo = dS.chi2Stereo_cut; % already delta-chi^2, but minimum not on grid point
%         Chi2KATRIN  =  chi2_KATRIN(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol); % absolute
%         Chi2Sum = Chi2Stereo + Chi2KATRIN-min(min(Chi2KATRIN));
%         DeltaChi2Combi_Shared = Chi2Sum - min(min(Chi2Sum));%
%         
%         % 2) KATRIN only
%         DeltaChi2Combi_KATRINonly = chi2_KATRIN-min(min(chi2_KATRIN(~dS.InterIdx)));
%         
%         % Combine 1)+2)
%         DeltaChi2Combi = zeros(1e3);
%         DeltaChi2Combi(dS.InterIdx) = DeltaChi2Combi_Shared;
%         DeltaChi2Combi(~dS.InterIdx) = DeltaChi2Combi_KATRINonly(~dS.InterIdx);
        
        Chi2Sum = dS.chi2Stereo+chi2_KATRIN;
        S.chi2_ref = min(min(Chi2Sum));
        S.chi2 = Chi2Sum;%DeltaChi2Combi;
        DeltaChi2Combi =  Chi2Sum-S.chi2_ref;
    else
        % sensitivity can be combined just as it is
%         Chi2Sum = dS.chi2Stereo +  chi2_KATRIN-min(min(chi2_KATRIN));       % absolute maps here
%         DeltaChi2Combi = Chi2Sum-min(min(Chi2Sum)); % minimum at null hypothesiis
    end
    %% get combined contour
    sin2T4_contour_Combi = zeros(1e3,1);
    mNu4Sq_contour_Combi = mNu4Sq(1,:);
    for i=1:1e3
        Idx = find(DeltaChi2Combi(:,i)>=dS.chi2Stereocrit(:,i),1,'first');
        if isempty(Idx)
            sin2T4_contour_Combi(i) = NaN;
        else
            sin2T4_contour_Combi(i) = sin2T4(Idx,i);
        end
    end

    save(savefile,...
        'sin2T4','mNu4Sq','DeltaChi2Combi',...
        'sin2T4_contour_STEREO','mNu4Sq_contour_STEREO',...
        'sin2T4_contour_KATRIN','mNu4Sq_contour_KATRIN',...
        'sin2T4_contour_Combi','mNu4Sq_contour_Combi',...
        'dS','chi2_KATRIN');
 
    if strcmp(DataType,'Real') 
        chi2_bf_Combi        = min(min(DeltaChi2Combi)); % 
        [row,col]            = find(DeltaChi2Combi==chi2_bf_Combi);
        sin2T4_bf_Combi      = sin2T4(row,col);
        mNu4Sq_bf_Combi      = mNu4Sq(row,col);
        sin2T4_bf_Combi      = sin2T4_bf_Combi;%(1:2);
        mNu4Sq_bf_Combi      = mNu4Sq_bf_Combi;%([1,3]);
        
        save(savefile,'sin2T4_bf_KATRIN','mNu4Sq_bf_KATRIN',...
            'sin2T4_bf_STEREO','mNu4Sq_bf_STEREO',...
            'sin2T4_bf_Combi','mNu4Sq_bf_Combi','-append');
    end
end

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);

%% Grid Plot
    DeltaChi2CombiPlt = DeltaChi2Combi;
    DeltaChi2CombiPlt(~dS.InterIdx & DeltaChi2CombiPlt>chi2inv(0.95,2))  = NaN;
    DeltaChi2CombiPlt(dS.InterIdx & DeltaChi2CombiPlt>dS.chi2Stereocrit) = NaN; 
    GetFigure;
  surf(sin2T4,mNu4Sq,DeltaChi2CombiPlt,'EdgeColor','none')
    view(2);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    grid off;
    c = colorbar;
    ylim([0.1 1600]);
    xlim([1e-03 0.5])
    PrettyFigureFormat('FontSize',22);
    c.Label.String = sprintf('\\Delta\\chi^2');
    c.FontSize = get(gca,'FontSize');
    ylim([0.1 38^2])
    xlim([3e-03 0.5])
    hold on;
    pC = plot3(sin2T4_contour_Combi,mNu4Sq_contour_Combi,5e4.*ones(1e3,1),'-','LineWidth',3,'Color',rgb('Black'));
    pK = plot3(sin2T4_contour_KATRIN,mNu4Sq_contour_KATRIN,5e4.*ones(numel(sin2T4_contour_KATRIN),1),'-.','LineWidth',2,'Color',rgb('DeepSkyBlue'));
    pS = plot3(sin2T4_contour_STEREO,mNu4Sq_contour_STEREO,5e4.*ones(1e3,1),':','LineWidth',2.5,'Color',rgb('Red'));
    
    ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
    xlabel(sprintf('|{\\itU_{e4}}|^2'));
   % 
    leg = legend([pK,pS,pC],...
        'KATRIN only','STEREO only','KATRIN + STEREO',...
        'Location','southwest');
    PrettyLegendFormat(leg);
    if strcmp(DataType,'Real')
        leg.Title.String = sprintf('Exclusion at 95%% C.L.');
    else
         leg.Title.String = sprintf('Sensitivity at 95%% C.L.');
    end
    leg.Title.FontWeight = 'normal';
    
    pltname1 = strrep(strrep(savefile,'results','plots'),'.mat','_GridPlot.png');
    print(gcf,pltname1,'-dpng','-r350');
    fprintf('save plot to %s \n',pltname1);
    
    return
%% Contour Plot
GetFigure;
pC = plot(sin2T4_contour_Combi,mNu4Sq_contour_Combi,'-','LineWidth',3,'Color',rgb('Black'));
hold on;
pK = plot(sin2T4_contour_KATRIN,mNu4Sq_contour_KATRIN,'--','LineWidth',2,'Color',rgb('DeepSkyBlue'));
pS = plot(sin2T4_contour_STEREO,mNu4Sq_contour_STEREO,':','LineWidth',2.5,'Color',rgb('Red'));
if strcmp(DataType,'Real')
    
    pCbf = plot(sin2T4_bf_Combi,mNu4Sq_bf_Combi,'s','LineWidth',3,'Color',pC.Color,'MarkerSize',8,'MarkerFaceColor',pC.Color);
    pKbf = plot(sin2T4_bf_KATRIN,mNu4Sq_bf_KATRIN,'x','LineWidth',2,'Color',pK.Color,'MarkerSize',8);
    pSbf = plot(sin2T4_bf_STEREO,mNu4Sq_bf_STEREO,'o','LineWidth',2,'Color',pS.Color,'MarkerSize',8,'MarkerFaceColor',pS.Color);
end
set(gca,'XScale','log');
set(gca,'YScale','log');
grid off;
ylim([0.1 1600]);
xlim([3e-03 0.5])
PrettyFigureFormat('FontSize',22);
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
xlabel(sprintf('|{\\itU_{e4}}|^2'));

leg = legend([pK,pS,pC],'KATRIN','STEREO','KATRIN + STEREO','Location','southwest');
if strcmp(DataType,'Real')
    leg.Title.String = sprintf('Exclusion at 95%% C.L.');
else
    leg.Title.String = sprintf('Sensitivity at 95%% C.L.');
end
leg.Title.FontWeight = 'normal';

PrettyLegendFormat(leg);
pltname2 = strrep(strrep(savefile,'results','plots'),'.mat','_ContourPlot.png');
print(gcf,pltname2,'-dpng','-r350');
fprintf('save plot to %s \n',pltname2);
%%
close all
GetFigure;
Chi2KATRIN  =  chi2_KATRIN(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol);
Chi2KATRIN = Chi2KATRIN -min(min(Chi2KATRIN));
Chi2KATRINPlt = Chi2KATRIN;
Chi2KATRINPlt(Chi2KATRINPlt>5.99) = NaN;
mK =surf(sin2T4(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol),mNu4Sq(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol),...
    Chi2KATRINPlt,'EdgeColor','none')
hold on;
[M,pContour] = contour3(sin2T4(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol),mNu4Sq(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol),...
    Chi2KATRINPlt,[0.1 0.5,0.65,0.7 2],'LineWidth',2,'LineStyle','-','Color',rgb('White'),...
                        'ShowText','off')
t =clabel(M,pContour,'FontSize',15,'Rotation',-45,'LineWidth',2,'Color',rgb('White'),'FontWeight','bold');
[~,cKATRIN] = contour3(sin2T4(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol),mNu4Sq(dS.Startrow:dS.Stoprow,dS.Startcol:dS.Stopcol),...
    Chi2KATRIN,[5.99,5.99],'LineWidth',2,'LineStyle','-','Color',rgb('Black'),...
                        'ShowText','off')
view(2)
set(gca,'XScale','log');
set(gca,'YScale','log');
grid off
c = colorbar;
PrettyFigureFormat('FontSize',22);
c.Label.String = sprintf('\\Delta\\chi^2');
c.FontSize = get(gca,'FontSize');
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
xlabel(sprintf('|{\\itU_{e4}}|^2'));
xlim([3e-03 0.5])
ylim([0.1 10])
if strcmp(DataType,'Real')
    leg = legend(cKATRIN,'KATRIN exclusion 95% C.L.','Location','southwest');
else
    leg = legend(cKATRIN,'KATRIN sensitivity 95% C.L.','Location','southwest');
end

PrettyLegendFormat(leg);
pltname2 = strrep(strrep(savefile,'results','plots'),'.mat','_GridKATRINonly.png');
print(gcf,pltname2,'-dpng','-r350');
fprintf('save plot to %s \n',pltname2);