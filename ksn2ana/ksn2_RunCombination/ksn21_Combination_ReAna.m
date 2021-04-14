% combine ksn1 and ksn2, nu-mass fixed
% simply add chi2-map2
% with KSN-1 Re-Ana chi2 map
%% settings that might change
chi2Name = 'chi2CMShape';
DataType = 'Real';
nGridSteps = 30;
range = 40;
BF = 'ON';
%% configure RunAnalysis object
if strcmp(chi2Name,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2Name,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar','E0 Norm Bkg',...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p1eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2Name,...
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
    'LoadGridArg',{'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'}};

%%
S = SterileAnalysis(SterileArg{:});
%% load ksn1
S.RunAnaObj.ModelObj.BKG_PtSlope = -2.2*1e-06;
S.RunAnaObj.DataSet = 'Knm1';
S.RunAnaObj.RunData.RunName = 'KNM1';
S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
S.RunAnaObj.AngularTFFlag ='ON';
S.RunAnaObj.SysBudget = 200;
S.RunAnaObj.FSDFlag = 'KNM2_0p1eV';
% stat and syst
S.RunAnaObj.NonPoissonScaleFactor = 1.064;
S.RunAnaObj.chi2 = 'chi2CMShape';
S.nGridSteps = 30;
S.LoadGridArg = {'CheckLarger','ON','ExtmNu4Sq','OFF','mNu4SqTestGrid',2};
S.LoadGridFile(S.LoadGridArg{:});
S.Interp1Grid('RecomputeFlag','ON');%,'Maxm4Sq',34.2^2);
p1tot = S.ContourPlot('BestFit',BF,'CL',95,'HoldOn','OFF','Color',rgb('FireBrick'),'LineStyle',':');

% load again with same binning as KSN-1 (restricted)
S.LoadGridFile(S.LoadGridArg{:}); 
S.Interp1Grid('RecomputeFlag','ON','Minm4Sq',1,'Maxm4Sq',38^2);
mNu4Sq_k1 = S.mNu4Sq;
sin2T4_k1 = S.sin2T4;
chi2_k1   = S.chi2;
chi2ref_k1= S.chi2_ref;
sum(sum(isnan(S.chi2)))
%% load ksn2
S.nGridSteps = 30;
S.RunAnaObj.ModelObj.BKG_PtSlope = 3*1e-06;
S.RunAnaObj.FSDFlag = 'KNM2_0p5eV';
S.RunAnaObj.DataSet = 'Knm2';
S.RunAnaObj.RunData.RunName = 'KNM2_Prompt';
S.RunAnaObj.ELossFlag  = 'KatrinT2A20';
S.RunAnaObj.AngularTFFlag ='ON';
S.RunAnaObj.SysBudget = 40;
S.RunAnaObj.NonPoissonScaleFactor = 1.112;

% stat and syst
S.RunAnaObj.chi2 = 'chi2CMShape';
S.LoadGridArg = {'mNu4SqTestGrid',5,'ExtmNu4Sq','ON'};
% load extended grid to get best fit
S.LoadGridFile(S.LoadGridArg{:},'Extsin2T4','ON');            
S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',35^2);
S.FindBestFit;
chi2_ref = S.chi2_ref;
mNuSq_bf = S.mNuSq_bf;
sin2T4_bf = S.sin2T4_bf;
mNu4Sq_bf = S.mNu4Sq_bf;
% load extended grid to get niceer plot
S.LoadGridFile(S.LoadGridArg{:},'Extsin2T4','OFF');  
S.Interp1Grid('RecomputeFlag','ON','Maxm4Sq',38^2);
S.chi2_ref  = chi2_ref;
S.mNuSq_bf  = mNuSq_bf ;
S.sin2T4_bf = sin2T4_bf;
S.mNu4Sq_bf = mNu4Sq_bf;
[p2tot,sinMin] = S.ContourPlot('BestFit',BF,'ReCalcBF','OFF','CL',95,'HoldOn','ON','Color',rgb('Orange'),'LineStyle','-.');


% mNu4Sq_k2 = S.mNu4Sq;
% sin2T4_k2 = S.sin2T4;
% chi2_k2   = S.chi2;
% chi2ref_k2= S.chi2_ref;
mNu4Sq_k2_contour = S.mNu4Sq_contour;
sin2T4_k2_contour = S.sin2T4_contour;

% load again with same binning as KSN-1 (restricted)
S.LoadGridFile(S.LoadGridArg{:}); 
S.Interp1Grid('RecomputeFlag','ON','Minm4Sq',1,'Maxm4Sq',38^2);
mNu4Sq_k2 = S.mNu4Sq;
sin2T4_k2 = S.sin2T4;
chi2_k2   = S.chi2;
chi2ref_k2= chi2_ref;

%% check if they have same binning
 if sum(sum(sin2T4_k2~=sin2T4_k1))>0 || sum(sum(mNu4Sq_k1~=mNu4Sq_k2))>0
     fprintf(2, 'KSN-1 and KSN-2 do not have  the same binning - return \n');
     return;
 end
%% combi 
% sum
chi2_sum = chi2_k1+chi2_k2;
S.chi2 = chi2_sum;
S.chi2_ref = chi2ref_k1+chi2ref_k2;
[p12tot,sinMin] = S.ContourPlot('BestFit','ON','CL',95,'HoldOn','ON','Color',rgb('DodgerBlue'),'LineStyle','-');

% find best fit;


leg = legend([p1tot,p2tot,p12tot],'KSN-1','KSN-2','KSN-1 and KSN-2 combined');
PrettyLegendFormat(leg);
% save
title([S.GetPlotTitle, sprintf(' , {\\itm}_\\nu^2 = 0 eV^2')],'FontWeight','normal')
xlim([2e-03,0.5]);
plotname = [extractBefore(S.DefPlotName,'Grid'),sprintf('KSN12_Combination_E0NormBkg_%s.png',chi2Name)];
print(gcf,plotname,'-dpng','-r300');



return;
%%  plot sum
DeltaChi2 = GetDeltaChi2(95,2);
chi2grid_k1 = chi2_k1;% S.chi2;
chi2grid_k1((chi2grid_k1-S.chi2_ref)>DeltaChi2) =  NaN;%DeltaChi2+chi2_ref;% NaN;
zlimMax = DeltaChi2;
GetFigure;

surf(sin2T4_k1,mNu4Sq_k1,chi2grid_k1-S.chi2_ref,'EdgeColor','interp','FaceColor','interp');
grid off
set(gca,'XScale','log')
set(gca,'YScale','log')
view([0 0 1])



%%
PlotGrid = 'OFF';
if strcmp(PlotGrid,'ON')
GetFigure;

% prepare
chi2grid_k1 = chi2_k1;
chi2grid_k1((chi2grid_k1-chi2ref_k1)>S.DeltaChi2) =  NaN;
chi2grid_k2 = chi2_k2;
chi2grid_k2((chi2grid_k2-chi2ref_k2)>S.DeltaChi2) =  NaN;
chi2grid_tot = chi2_sum;
chi2ref_tot = chi2ref_k1+chi2ref_k2;
chi2grid_tot((chi2grid_tot-chi2ref_tot)>S.DeltaChi2) =  NaN;

% plot
%surf(sin2T4_k1,mNu4Sq_k1,chi2grid_k1-chi2ref_k1,'EdgeColor','interp','FaceColor','interp')
%hold on;
%surf(sin2T4_k2,mNu4Sq_k2,chi2grid_k2-chi2ref_k2,'EdgeColor','interp','FaceColor','interp')
surf(sin2T4_k2,mNu4Sq_k2,chi2grid_tot-chi2ref_tot,'EdgeColor','interp','FaceColor','interp')


PrettyFigureFormat('FontSize',20);
zlimMax = S.DeltaChi2;
zlim([0 zlimMax])
set(gca,'XScale','log')
set(gca,'YScale','log')
c =colorbar;
c.Label.String = sprintf('\\Delta\\chi^2');
c.Label.FontSize = get(gca,'FontSize')+2;
c.Limits=[0 zlimMax];
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
zlabel(sprintf('\\Delta\\chi^2'))
grid off
view([0 0 1]);


xlim([1e-03,0.5]);

end

%%
%% load ksn1
S.nGridSteps = 50;

S.RunAnaObj.DataSet = 'Knm1';
S.RunAnaObj.RunData.RunName = 'KNM1';
S.RunAnaObj.ELossFlag  = 'KatrinT2';
S.RunAnaObj.AngularTFFlag ='OFF';
S.RunAnaObj.SysBudget = 24;

% stat and syst
S.RunAnaObj.NonPoissonScaleFactor = 1.064;
S.RunAnaObj.chi2 = 'chi2CMShape';
S.LoadGridFile('CheckLarger','ON');
S.Interp1Grid('RecomputeFlag','ON');

hold on;
%%
p2 =plot3(sin2T4_k2_contour, mNu4Sq_k2_contour,S.DeltaChi2.*ones(size(mNu4Sq_k2_contour)));
