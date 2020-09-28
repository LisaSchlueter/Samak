%% load or calc
range = 40;
AuxLines = 'OFF';
ShowResults = 'OFF';
SavePlot = 'ON';
AnaMode = 'Uniform';%MultiRing';
SysBudget = 38;
DataType = 'Real';
if ~strcmp(AnaMode,'MultiRing')
    AnaStr= '';
else
    AnaStr = '_MultiRing';
end
%%
savedir = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/'];

savenameStat = sprintf('%sknm2ub1_Chi2Curve_%s_%.0feV_%s_chi2Stat.mat',...
    savedir,DataType,range,strrep(freePar,' ',''));
savenameCM = sprintf('%sknm2ub1_Chi2Curve_%s_%.0feV_%s_chi2CMShape_SysBudget%.0f.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),SysBudget);
try
dStat = importdata(savenameStat);
dCM = importdata(savenameCM);
catch
    fprintf('files not available \n')
    fprintf('%s \n',savenameCM)
    return
end

%% prepare vars: 
%stat
ScanResultStat = dStat.ScanResult;
FitResultStat = dStat.FitResult;
nFitMax = size(ScanResultStat.chi2min,1);
mNuSqStat  = reshape(ScanResultStat.ParScan,[nFitMax*2,1]);
mNuSqStat = mNuSqStat(~isnan(mNuSqStat));
[mNuSqStat,sortI] = sort(mNuSqStat);
mychi2minStat = reshape(ScanResultStat.chi2min,[nFitMax*2,1]);
mychi2minStat = mychi2minStat(~isnan(mychi2minStat));
mychi2minStat = mychi2minStat(sortI);
% stat + syst
ScanResultCM = dCM.ScanResult;
FitResultCM = dCM.FitResult;
nFitMax = size(ScanResultCM.chi2min,1);
mNuSqCM  = reshape(ScanResultCM.ParScan,[nFitMax*2,1]);
mNuSqCM  = mNuSqCM(~isnan(mNuSqCM));
[mNuSqCM,sortI] = sort(mNuSqCM);
mychi2minCM = reshape(ScanResultCM.chi2min,[nFitMax*2,1]);
mychi2minCM = mychi2minCM(~isnan(mychi2minCM));
mychi2minCM = mychi2minCM(sortI);
%% plot

 f4 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
pCM = plot(mNuSqCM,mychi2minCM,'-','LineWidth',3,'Color',rgb('DodgerBlue'));
hold on;
pStat = plot(mNuSqStat,mychi2minStat,'-.','LineWidth',3,'Color',rgb('SkyBlue'));


if strcmp(AuxLines,'ON')
p1 =plot(FitResultStat.par(1)+ScanResultStat.AsymErr(1).*ones(100,1),linspace(0,1,100),...
    ':','LineWidth',3,'Color',rgb('GoldenRod'));
p2 =plot(FitResultStat.par(1)+ScanResultStat.AsymErr(2).*ones(100,1),linspace(0,1,100),...
    ':','LineWidth',3,'Color',rgb('GoldenRod'));
p3 = plot(linspace(ScanResultStat.AsymErr(2)+FitResultStat.par(1),...
    ScanResultStat.AsymErr(1)+FitResultStat.par(1),100),ones(1,100),...
    ':','LineWidth',3,'Color',rgb('GoldenRod'));
end

xstr = sprintf('{{\\itm}_\\nu}^2');
xUnit = sprintf('eV^2');
PrettyFigureFormat('FontSize',24);

xlabel(sprintf('%s (%s)',xstr,xUnit),'Interpreter','tex');
ylabel(sprintf('\\chi^2 (%.0f dof)',ScanResultStat.dof(1,1) - 1 ));

if abs(FitResultStat.par(1))>0.05
    parStr = sprintf('%.3f',FitResultStat.par(1));
else
    parStr = sprintf('%.3f',FitResultStat.par(1));
end
if strcmp(ShowResults,'ON')
leg = legend([pStat,pCM],...
   sprintf('Stat. only:        %s = %.3f (%.3f +%.3f) %s',xstr,FitResultStat.par(1), ScanResultStat.AsymErr(2),ScanResultStat.AsymErr(1),xUnit),...
   sprintf('Stat. and syst: %s = %.3f (%.3f +%.3f) %s',xstr,FitResultCM.par(1), ScanResultCM.AsymErr(2),ScanResultCM.AsymErr(1),xUnit));
else
    leg = legend([pStat,pCM],'Stat. only','Stat. and syst.');
end
leg.EdgeColor = rgb('Silver');
leg.Location = 'northwest';
xlim([-0.8,0.8]);
switch DataType
    case 'Twin'
t = title('KNM2 twins','FontWeight','normal','FontSize',get(gca,'FontSize'));
    case 'Real'
        t = title(sprintf('KNM2 data'),'FontWeight','normal','FontSize',get(gca,'FontSize'));
end
%% save
if strcmp(SavePlot,'ON')
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = sprintf('%sknm2FS2_Chi2CurveStatSyst%s.pdf',plotdir,AnaStr);
fprintf('save plot to %s \n',plotname);
export_fig(gcf,plotname);
end

%% 
fprintf('--------%s----------------------------------\n',AnaMode)
fprintf('stat. only mnu^2 = %.3f (-%.3f + %.3f) eV^2 \n',FitResultStat.par(1),-FitResultStat.errNeg(1),FitResultStat.errPos(1));
fprintf('syst. only mnu^2 = %.3f (-%.3f + %.3f) eV^2 \n',FitResultCM.par(1),sqrt(FitResultCM.errNeg(1)^2-FitResultStat.errNeg(1)^2),sqrt(FitResultCM.errPos(1)^2-FitResultStat.errPos(1)^2));
fprintf('total      mnu^2 = %.3f (-%.3f + %.3f) eV^2 \n',FitResultCM.par(1),-FitResultCM.errNeg(1),FitResultCM.errPos(1));
