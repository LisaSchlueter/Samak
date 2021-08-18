% ksn2 calculate chi2-prifle for nu-mass
% perform grid searches for differed nu-masses
% exclude large mixings
%% settings that might change
DataType = 'Twin';
switch DataType
    case 'Twin'
        FixmNuSq_all = -1.1:0.1:2.5;
    case 'Real'
        FixmNuSq_all = sort(round([-3:0.1:2.5],2));
end

RecomputeFlag = 'OFF';
Maxm4Sq    = 36^2; % interpolation
chi2min = zeros(numel(FixmNuSq_all),1);
sin2T4_bf =  zeros(numel(FixmNuSq_all),1);
mNu4Sq_bf =  zeros(numel(FixmNuSq_all),1);

savedir       = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefileCombi = sprintf('%sksn2_%s_NuMassSensitivityGridSearch_ExclLargeUe4.mat',...
    savedir,DataType);

if exist(savefileCombi,'file') && strcmp(RecomputeFlag,'OFF')
    load(savefileCombi);
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
        'fixPar','E0 Norm Bkg',...%free par
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
    switch DataType
        case 'Real'
            LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON','CheckSmallerN','ON'};
            nGridSteps = 40;
        case 'Twin'
            LoadGridArg = {'mNu4SqTestGrid',2,'ExtmNu4Sq','ON'};
            nGridSteps = 30;
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
      
    ExclUe4 = sort([0.5,0.4:-0.1:0.1,0.01,0.02,0.03,0.04,0.05,0.07,0.005,0.002],'descend');%  
    
    chi2min_Combi   = zeros(numel(FixmNuSq_all),numel(ExclUe4));
    sin2T4_Combi_bf =  zeros(numel(FixmNuSq_all),numel(ExclUe4));
    mNu4Sq_Combi_bf =  zeros(numel(FixmNuSq_all),numel(ExclUe4));
    
    %%
    for i=1:numel(FixmNuSq_all)
        FixmNuSq = FixmNuSq_all(i);
        A.ModelObj.mnuSq_i = FixmNuSq;
        S.LoadGridFile(S.LoadGridArg{:},'FixmNuSq',FixmNuSq);
        S.Interp1Grid('MaxM4Sq',Maxm4Sq);
        
        if S.chi2_bf<S.chi2_Null
            chi2min(i) = S.chi2_bf;
            mNu4Sq_bf(i) = S.mNu4Sq_bf;
            sin2T4_bf(i) = S.sin2T4_bf;
        else
            chi2min(i) = S.chi2_Null;
            mNu4Sq_bf(i) = 0;
            sin2T4_bf(i) = 0;
        end
        
        chi2_KATRIN = S.chi2;
        for j=1:numel(ExclUe4)
            chi2Sum = chi2_KATRIN;
            chi2Sum(S.sin2T4>=ExclUe4(j)) = 1e3;
            
            S.chi2_ref = min(min(chi2Sum));
            S.chi2 = chi2Sum;
%             S.GridPlot('BestFit','ON','Contour','ON','CL',99.99,...
%                 'SavePlot','OFF',...
%                 'ExtraStr',sprintf('STEREOCombi_mNuSq%.2geV2',FixmNuSq));
%             close;
              S.FindBestFit;
            chi2min_Combi(i,j) = S.chi2_bf;
            mNu4Sq_Combi_bf(i,j) = S.mNu4Sq_bf;
            sin2T4_Combi_bf(i,j) = S.sin2T4_bf;
        end
        
    end
    save(savefileCombi,'chi2min','mNu4Sq_bf','sin2T4_bf','FixmNuSq_all',...
        'chi2min_Combi','mNu4Sq_Combi_bf','sin2T4_Combi_bf','ExclUe4')
    
end
 
%% get 1 sigma err
DeltaChi2 = chi2min_Combi-min(chi2min_Combi);
[MinIdx,~] = find(DeltaChi2==0);

ErrLow     = zeros(numel(ExclUe4),1);
ErrUp      = zeros(numel(ExclUe4),1);
ErrLow2    = zeros(numel(ExclUe4),1);
ErrUp2     = zeros(numel(ExclUe4),1);

%mNuSq_bf  = zeros(numel(ExclUe4),1);
%if strcmp(DataType,'Real')
    mNu_tmp  = linspace(min(FixmNuSq_all),max(FixmNuSq_all),1e3);
    chi2_tmp = interp1(FixmNuSq_all,chi2min_Combi,mNu_tmp,'spline');
    [Idx,row]= find(chi2_tmp==min(chi2_tmp));
    mNuSq_bf = mNu_tmp(Idx);
    deltachi2_tmp = chi2_tmp-min(chi2_tmp);
%end

for i=1:numel(ExclUe4)
ErrLow(i) = mNuSq_bf(i)-interp1(deltachi2_tmp(1:Idx(i),i),mNu_tmp(1:Idx(i)),1,'spline');
ErrUp(i)  = interp1(deltachi2_tmp(Idx(i):end,i),mNu_tmp(Idx(i):end),1,'spline')-mNuSq_bf(i);
% 2 sigma
ErrLow2(i) = mNuSq_bf(i)-interp1(deltachi2_tmp(1:Idx(i),i),mNu_tmp(1:Idx(i)),4,'spline');
ErrUp2(i)  = interp1(deltachi2_tmp(Idx(i):end,i),mNu_tmp(Idx(i):end),4,'spline')-mNuSq_bf(i);
end


%%
PltIdx = find(sum(ExclUe4==[0.5,0.4,0.2,0.1,0.01,0.02,0.04,0.05,0.002]'));

pHandle = cell(numel(PltIdx),1);
legStr = cell(numel(PltIdx),1);
myColors = parula(numel(PltIdx));
LineStyles = {'-.',':','--','-','-.',':','--','-','-.',':','--','-','-.',':','--'};
figure('Units','normalized','Position',[0.1,0.1,0.8,0.5]);
plot(FixmNuSq_all,ones(numel(FixmNuSq_all),1),'k-','LineWidth',2);
hold on;
plot(FixmNuSq_all,4.*ones(numel(FixmNuSq_all),1),'k-','LineWidth',2);
for i=1:numel(PltIdx)
pHandle{i} = plot(FixmNuSq_all,DeltaChi2(:,PltIdx(i)),'Color',myColors(i,:),'LineStyle',LineStyles{i},'LineWidth',3);
if abs(ErrLow(PltIdx(i)))>100
    if strcmp(DataType,'Twin')
        legStr{i} = sprintf('|{\\itU}_{e4}|^2 < %.3f eV^2 , \\sigma({\\itm}_\\nu^2) = - \\infty + %.2f eV^2',...
            ExclUe4(PltIdx(i)),abs(ErrUp(PltIdx(i))));
    else
        legStr{i} = sprintf('|{\\itU}_{e4}|^2 < %.3f eV^2 ,  {\\itm}_\\nu^2 = %.2f (-\\infty    + %.2f) eV^2',...
            ExclUe4(PltIdx(i)),mNuSq_bf(PltIdx(i)),abs(ErrUp(PltIdx(i))));
    end
else
    if strcmp(DataType,'Twin')
        legStr{i} = sprintf('|{\\itU}_{e4}|^2 < %.3f eV^2 ,  \\sigma({\\itm}_\\nu^2) = -%.2f + %.2f eV^2',...
            ExclUe4(PltIdx(i)),abs(ErrLow(PltIdx(i))),abs(ErrUp(PltIdx(i))));
    else
        legStr{i} = sprintf('|{\\itU}_{e4}|^2 < %.3f eV^2  , {\\itm}_\\nu^2 = %.2f (-%.2f + %.2f) eV^2',...
            ExclUe4(PltIdx(i)),mNuSq_bf(PltIdx(i)) ,  abs(ErrLow(PltIdx(i))),abs(ErrUp(PltIdx(i))));
    end
end

end
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('\\Delta\\chi^2'));
PrettyFigureFormat('FontSize',24);

%% plot knm2 chi2 profile

PltKnm2 = 'ON';
if strcmp(PltKnm2,'ON')
    dN = importdata(sprintf('%stritium-data/fit/Knm2/Chi2Profile/Uniform/Chi2Profile_%s_UniformScan_mNu_Knm2_UniformFPD_chi2CMShape_SysBudget40_NP1.112_FitParE0BkgNorm_nFit50_KNM2_0p1eV_min-2.5_max2.5.mat',getenv('SamakPath'),DataType));
    pKnm2 = plot(dN.ScanResults.ParScan(:,1),dN.ScanResults.chi2min(:,1)-min(min(dN.ScanResults.chi2min)),'-','Color',rgb('FireBrick'),'LineWidth',3);
    plot(dN.ScanResults.ParScan(1:end,2),dN.ScanResults.chi2min(1:end,2)-min(min(dN.ScanResults.chi2min)),'-','Color',pKnm2.Color,'LineWidth',pKnm2.LineWidth);
    if strcmp(DataType,'Twin')
         knm2leg = sprintf('No sterile neutrinos, \\sigma({\\itm}_\\nu^2) = -%.2f + %.2f eV^2',...
        dN.ScanResults.BestFit.errNeg,dN.ScanResults.BestFit.errPos);
    else
    knm2leg = sprintf('No sterile neutrinos, {\\itm}_\\nu^2 = %.2f (-%.2f + %.2f) eV^2',...
        dN.ScanResults.BestFit.par(1),dN.ScanResults.BestFit.errNeg,dN.ScanResults.BestFit.errPos);
    end
    leg = legend([pHandle{:},pKnm2],{legStr{:},knm2leg});
else
    leg = legend([pHandle{:}],legStr);
end
%%

PrettyLegendFormat(leg);

leg.NumColumns = 1;
leg.Location = 'eastoutside';
leg.FontSize =  get(gca,'FontSize')-4;
if strcmp(DataType,'Twin')
    ylim([0 7]);
    xlim([min(FixmNuSq_all),1.5])
    leg.Title.String = sprintf('MC twin + constraint on |{\\itU}_{e4}|^2');
else
  
    ylim([0 7]);
    xlim([-1.5,2.5])
    leg.Title.String = sprintf('Data + constraint on |{\\itU}_{e4}|^2');

end
leg.Title.FontWeight = 'normal';
leg.Title.FontSize = get(gca,'FontSize');
%t.FontWeight = 'normal'; t.FontSize = get(gca,'FontSize');

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = strrep(strrep(savefileCombi,'results','plots'),'.mat','.pdf');
%print(gcf,pltname,'-dpng','-r350');
export_fig(pltname);
fprintf('save plot to %s \n',pltname);
