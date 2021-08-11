
% combine ksn1 and ksn2, nu-mass free
% plot with ksn1, ksn2, ksn1+2
chi2          = 'chi2CMShape';
DataType      = 'Real';
nGridStepsCommon  = 30;
freePar       = 'mNu E0 Bkg Norm';
range         = 40;
extraStr = '_ReAna';
FixBfK2 = 'OFF';
N4 = 'OFF';
ExtLeg = 'ON';
savedir = [getenv('SamakPath'),'/SterileAnalysis/GridSearchFiles/Combi/',DataType,'/'];
MakeDir(savedir)
savename = sprintf('%sKSN12Combi%s_GridSearch_%s_%s_Uniform_%s_%.0fnGrid.mat',...
    savedir,extraStr,DataType,strrep(freePar,' ',''),chi2,nGridStepsCommon);%_ReAna

% load combi file
if exist(savename,'file') 
    d = importdata(savename);
    fprintf('load results from file %s \n',savename)
else
     fprintf('file doesnt exist %s \n',savename)
    return
end

savename2 = sprintf('%sKSN12Combi%s_RunAnaObj_%s.mat',savedir,extraStr,DataType);
if exist(savename2,'file')
    A = importdata(savename2);
else
    % int obejct for interpolatio and plot
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','mNu E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p1eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'NonPoissonScaleFactor',1.112,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    save(savename2,'A');
end

A.fixPar = ConvertFixPar('freePar','mNu E0 Norm Bkg','nPar',A.nPar,'nPixel',1);
A.chi2 = chi2;

%% configure Sterile analysis object
SterileArg = {'RunAnaObj',A,... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',50,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'RandMC','OFF',...
    'range',range};
S = SterileAnalysis(SterileArg{:});
S.NullHypothesis = 'OFF';
%% ksn1
KSN1config = 'New';
S.RunAnaObj.DataSet = 'Knm1';
S.RunAnaObj.RunData.RunName = 'KNM1';
S.RunAnaObj.NonPoissonScaleFactor = 1.064;
S.RunAnaObj.chi2 = 'chi2CMShape';
switch KSN1config
    case 'Old'
        % OLD settings (PRL publication)
        S.nGridSteps = 50;
        S.RunAnaObj.ELossFlag  = 'KatrinT2';
        S.RunAnaObj.AngularTFFlag ='OFF';
        S.RunAnaObj.SysBudget = 24;
        S.RunAnaObj.FSDFlag = 'Sibille0p5eV';
        S.RunAnaObj.ModelObj.BKG_PtSlope = 0;
        S.LoadGridArg = '';
        S.InterpMode = 'Mix';
        % load map
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',19^2);
    case 'New'
        % Re-Analysis
        S.nGridSteps = 50;
        S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
        S.RunAnaObj.AngularTFFlag ='ON';
        S.RunAnaObj.SysBudget = 200;
        S.RunAnaObj.FSDFlag = 'KNM2_0p1eV';
        S.RunAnaObj.ModelObj.BKG_PtSlope = -2.2*1e-06;
        S.LoadGridArg = {'mNu4SqTestGrid',2};
        S.InterpMode = 'Mix';
        % load map
        S.LoadGridFile(S.LoadGridArg{:});
        S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',36^2);
end

if strcmp(DataType,'Real')
    BF = 'ON';
else
    BF = 'OFF';
end

% store
mNu4Sq_k1 = S.mNu4Sq;
sin2T4_k1 = S.sin2T4;
chi2_k1   = S.chi2;
chi2ref_k1= S.chi2_ref;
sum(sum(isnan(S.chi2)))

% plot
[p1tot,~,pbf1] = S.ContourPlot('BestFit',BF,'CL',95,'HoldOn','OFF','Color',rgb('FireBrick'),'LineStyle',':','MarkerStyle','*');
if strcmp(DataType,'Real')
    mNu4Sq_1_contour = S.mNu4Sq_contour;
    sin2T4_k1_contour = S.sin2T4_contour;
else
    mNu4Sq_1_contour = S.mNu4Sq_contour(:,3);
    sin2T4_k1_contour = S.sin2T4_contour(:,3);
end

dof1 = S.dof;

% find significance 
[DeltaChi2_1, SignificanceBF_1,SignificanceSigma_1] = S.CompareBestFitNull;
mNu4Sqbf_k1 = S.mNu4Sq_bf;
sin2T4bf_k1 = S.sin2T4_bf;
%% load ksn2
S.InterpMode = 'spline';
if strcmp(DataType,'Real')
S.nGridSteps = 50;
else
    S.nGridSteps = 40;
end
S.RunAnaObj.FSDFlag = 'KNM2_0p5eV';
S.RunAnaObj.DataSet = 'Knm2';
S.RunAnaObj.RunData.RunName = 'KNM2_Prompt';
S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
S.RunAnaObj.AngularTFFlag ='ON';
S.RunAnaObj.SysBudget = 40;
S.RunAnaObj.NonPoissonScaleFactor = 1.112;
S.RunAnaObj.ModelObj.BKG_PtSlope = 3.*1e-06;
% stat and syst
S.RunAnaObj.chi2 = 'chi2CMShape';
S.LoadGridFile('CheckLarger','OFF','mNu4SqTestGrid',5,'ExtmNu4Sq','OFF');            
S.Interp1Grid('RecomputeFlag','ON');
mNu4Sq_k2 = S.mNu4Sq;
sin2T4_k2 = S.sin2T4;
chi2_k2   = S.chi2;
chi2ref_k2= S.chi2_ref;
mNu4Sq_k2_contour = S.mNu4Sq_contour;
sin2T4_k2_contour = S.sin2T4_contour;
%% get best fit from fine grid
savedir = sprintf('%sksn2ana/ksn2_BestFit/results/',getenv('SamakPath'));
savefile2 = sprintf('%sksn2_ImpBestFit_%s_%s_%s_nGridStepsFull%.0f_nGridStepsImp%.0f.mat',...
    savedir,DataType,'mNuE0NormBkg',chi2,50,40);
if exist(savefile2,'file') 
    d2 = importdata(savefile2);
    S.mNu4Sq_bf = d2.mNu4Sq_bf_inter;
    S.sin2T4_bf = d2.sin2T4_bf_inter;
    S.chi2_bf = d2.chi2_bf_inter;
else
    S.FindBestFit;
end
mNu4Sqbf_k2 = S.mNu4Sq_bf;
sin2T4bf_k2 = S.sin2T4_bf;

[p2tot,sinMin,pbf2] = S.ContourPlot('BestFit',BF,'CL',95,'HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.','ReCalcBF','OFF');
dof2 = S.dof;
% find significance
[DeltaChi2_2, SignificanceBF_2,SignificanceSigma_2] = S.CompareBestFitNull;
%% KSN1+2
S.sin2T4 = d.sin2T4;
S.mNu4Sq = d.mnu4Sq;
S.chi2 = d.chi2;
S.chi2_Null = d.FitResult_Null.chi2min;
S.mNuSq = cell2mat(cellfun(@(x) x.par(1),d.FitResults,'UniformOutput',0));
S.E0 = cell2mat(cellfun(@(x) x.par(2),d.FitResults,'UniformOutput',0));
dof12 = 46;
S.dof = dof12;
if strcmp(DataType,'Twin') % 2 fits didn't converge....
    S.chi2(21,11) = NaN;
    S.chi2(25,11) = NaN;
end
S.InterpMode = 'spline';
S.Interp1Grid('Maxm4Sq',40^2);%,'Minm4Sq',40);
  
if strcmp(FixBfK2,'ON')
    % fix location of best-fit to ksn-2 only
%     IdxSin = find(abs(sin2T4bf_k2-S.sin2T4(:,1))==min(abs(sin2T4bf_k2-S.sin2T4(:,1))));
%     IdxM4 = find(abs(mNu4Sqbf_k2-S.mNu4Sq(1,:))==min(abs(mNu4Sqbf_k2-S.mNu4Sq(1,:))));
    S.chi2_ref = chi2ref_k1+chi2ref_k2+2;%S.chi2(IdxSin,IdxM4);
    BF = 'OFF';
elseif strcmp(S.NullHypothesis,'ON')
    S.chi2_ref = d.FitResult_Null.chi2min;
else
    S.FindBestFit;
   % S.chi2_ref = d.chi2_ref;
end
%

[p12tot,~,pbf12] = S.ContourPlot('HoldOn','ON','BestFit',BF,...
                                'SavePlot','OFF','ReCalcBf','OFF','MarkerStyle','o');

mNu4Sq_k12_contour = S.mNu4Sq_contour;
sin2T4_k12_contour = S.sin2T4_contour;

mNu4Sqbf_k12 = S.mNu4Sq_bf;
sin2T4bf_k12 = S.sin2T4_bf;

% find significance
[DeltaChi2_tot, SignificanceBF_tot,SignificanceSigma_tot] = S.CompareBestFitNull;

% convert cl to sigma
SigmaBF1 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_1,'nPar',2);
SigmaBF2 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_2,'nPar',2);
SigmaBF12 = ConvertCLStd('Mode','CL2Sigma','CL',SignificanceBF_tot,'nPar',2);
%%

legStr = sprintf('%s: {\\itm}_\\nu^2 free at %.0f%% C.L.',S.GetPlotTitle('Mode','data'),S.ConfLevel);
if strcmp(ExtLeg,'ON')
    leg = legend([p1tot,p2tot,p12tot,pbf1,pbf2,pbf12],...
        'KSN-1','KSN-2','KSN-1 and KSN-2',...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sqbf_k1,sin2T4bf_k1,SigmaBF1),...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sqbf_k2,sin2T4bf_k2,SigmaBF2),...
        sprintf('{\\itm}_4^2 = %.1f eV^2 , |{\\itU}_{e4}|^2 = %.3f (%.1f \\sigma)',mNu4Sqbf_k12,sin2T4bf_k12,SigmaBF12));
    leg.NumColumns = 2;
    title(legStr,'FontWeight','normal','FontSize',get(gca,'FontSize'));
else
    leg = legend([p1tot,p2tot,p12tot],'KSN-1','KSN-2','KSN-1 and KSN-2 combined');
    leg.Title.String = legStr;
    leg.Title.FontWeight = 'normal';
    title('');
end
PrettyLegendFormat(leg);

xlim([2e-03,0.5]);


%%
if strcmp(N4,'ON')
    extraStr = '_N4';
    savedirOther = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];
    filenameN4 = sprintf('%scoord_Neutrino4_123sigma.mat',savedirOther);
    dN4 = importdata(filenameN4);
    % convet 2 tritium par space
    sin2T4_N4 = sin(asin(sqrt(dN4.SinSquare2Theta_X_2sigma))./2).^2;
    pN4 = plot(sin2T4_N4,dN4.DmSquare41_Y_2sigma,'k-','LineWidth',2);
    leg = legend([p1tot,p2tot,p12tot,pN4],'KSN-1','KSN-2','KSN-1 and KSN-2 combined',sprintf('Neutrino-4 (2\\sigma)'));
    PrettyLegendFormat(leg);
     leg.NumColumns = 1;
else
    extraStr = '';
    
end

 if strcmp(ExtLeg,'ON')
    ylim([0.5 42^2]);
else
    ylim([3 40^2]);
end
%%
% save plot
plotname = [extractBefore(S.DefPlotName,'Grid'),sprintf('KSN12_%s_Combination_mNuE0NormBkg_%s%s.png',...
    S.GetPlotTitle('Mode','data'),chi2,extraStr)];
if strcmp(ExtLeg,'ON')
    plotname = strrep(plotname,'.png','_BF.png');
end
print(gcf,plotname,'-dpng','-r300');
fprintf('save plot to %s \n',plotname);

export_fig(gcf,strrep(plotname,'.png','.pdf'));

%% save combi file
savedir = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));
MakeDir(savedir);
savefile = sprintf('%sksn21_Combi_freemNuSq_ReAna_%s.mat',savedir,DataType);
save(savefile,'chi2ref_k1','chi2ref_k2','DeltaChi2_1','DeltaChi2_2','DeltaChi2_tot',...
    'mNu4Sq_1_contour','mNu4Sq_k2_contour','sin2T4_k1_contour','sin2T4_k2_contour',...
    'mNu4Sq_k12_contour','sin2T4_k12_contour','d',...
    'mNu4Sqbf_k1','mNu4Sqbf_k2','mNu4Sqbf_k12',...
    'sin2T4bf_k1','sin2T4bf_k2','sin2T4bf_k12',...
    'dof1','dof2','dof12');
