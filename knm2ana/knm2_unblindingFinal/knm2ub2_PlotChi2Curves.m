%% load or calc
range = 40;
AuxLines = 'ON';
ShowResults = 'ON';
SavePlot = 'ON';
AnaFlag = 'StackPixel';%StackPixel';

if strcmp(AnaFlag,'Ring')
    SysBudget = 39;
else
    SysBudget = 38;
end
DataType = 'Real';
freePar = 'mNu E0 Bkg Norm';
%%
savedir = [getenv('SamakPath'),'knm2ana/knm2_unblindingFinal/results/Chi2Curve/'];

savenameStat = sprintf('%sknm2ub2_Chi2Curve_%s_%.0feV_%s_chi2Stat_%s.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),AnaFlag);
savenameCM = sprintf('%sknm2ub2_Chi2Curve_%s_%.0feV_%s_chi2CMShape_%s_SysBudget%.0f.mat',...
    savedir,DataType,range,strrep(freePar,' ',''),AnaFlag,SysBudget);
try
dStat = importdata(savenameStat);
fprintf('load %s  \n',savenameStat);
dCM = importdata(savenameCM);
fprintf('load %s  \n',savenameCM);
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
 if strcmp(AuxLines,'ON')
% p1 =plot(FitResultStat.par(1)+ScanResultStat.AsymErr(1).*ones(100,1),linspace(0,1,100),...
%     ':','LineWidth',3,'Color',rgb('GoldenRod'));
% p2 =plot(FitResultStat.par(1)+ScanResultStat.AsymErr(2).*ones(100,1),linspace(0,1,100),...
%     ':','LineWidth',3,'Color',rgb('GoldenRod'));
% p3 = plot(linspace(ScanResultStat.AsymErr(2)+FitResultStat.par(1),...
%     ScanResultStat.AsymErr(1)+FitResultStat.par(1),100),ones(1,100),...
%     ':','LineWidth',3,'Color',rgb('GoldenRod'));


ConfICM = area(dCM.FitResult.par(1)+[FitResultCM.errNeg(1),FitResultCM.errPos(1)],[200,200],...
    'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.2,'EdgeColor','none');
hold on;
ConfIStat = area(dStat.FitResult.par(1)+[FitResultStat.errNeg(1),FitResultStat.errPos(1)],[200,200],...
    'FaceColor',rgb('Orange'),'FaceAlpha',0.25,'EdgeColor','none');

plot(dStat.FitResult.par(1).*ones(100,1),linspace(0,2e2,1e2),'--','LineWidth',2,'Color',rgb('Orange'))
plot(dCM.FitResult.par(1).*ones(100,1),linspace(0,2e2,1e2),'--','LineWidth',2,'Color',rgb('DodgerBlue'))
 end

pCM = plot(mNuSqCM,mychi2minCM,'-','LineWidth',3,'Color',rgb('DodgerBlue'));
hold on;
pStat = plot(mNuSqStat,mychi2minStat,'-','LineWidth',3,'Color',rgb('Orange'));


xstr = sprintf('{{\\itm}_\\nu}^{2}');
xUnit = sprintf('eV^{ 2}');
PrettyFigureFormat('FontSize',22);

xlabel(sprintf('%s (%s)',xstr,xUnit),'Interpreter','tex');
ylabel(sprintf('\\chi^2 (%.0f dof)',ScanResultStat.dof(1,1) - 1 ));

if abs(FitResultStat.par(1))>0.05
    parStr = sprintf('%.3f',FitResultStat.par(1));
else
    parStr = sprintf('%.3f',FitResultStat.par(1));
end
if strcmp(ShowResults,'ON')
leg = legend([pStat,pCM],...
   sprintf('Stat. only:        %s = %.2f (+%.2f %.2f) %s ',xstr,FitResultStat.par(1),FitResultStat.errPos(1),FitResultStat.errNeg(1),xUnit),...
   sprintf('Stat. and syst: %s = %.2f (+%.2f %.2f) %s ',xstr,FitResultCM.par(1),FitResultCM.errPos(1),FitResultCM.errNeg(1),xUnit));
else
    leg = legend([pStat,pCM],'Stat. only','Stat. and syst.');
end
leg.EdgeColor = rgb('Silver');
leg.LineWidth  = 1;
leg.Location = 'north';
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.9]));
if strcmp(AnaFlag,'StackPixel')
    AnaStr = 'Uniform';
else
    AnaStr = 'Multi-Ring (4)';
end

switch DataType
    case 'Twin'
        xlim([-1,1]);
        ylim([0 11]);
        t = title(sprintf('KNM2 twins - %s',AnaStr),'FontWeight','normal','FontSize',get(gca,'FontSize'));
    case 'Real'
        if strcmp(AnaFlag,'Ring')
            xlim([-1.15,0.7]);
            ylim([dCM.FitResult.chi2min-1 1+max(max(dStat.ScanResult.chi2min))])
        else
            xlim([-0.6,1]);
            ylim([dCM.FitResult.chi2min-1 1+max(max(dCM.ScanResult.chi2min))])
        end
        t = title(sprintf('KNM2 data - %s',AnaStr),'FontWeight','normal','FontSize',get(gca,'FontSize'));
end
%% save
if strcmp(SavePlot,'ON')
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = sprintf('%sknm2ub2_Chi2CurveStatSyst_%s_%s.pdf',plotdir,strrep(freePar,' ',''),AnaFlag);
fprintf('save plot to %s \n',plotname);
export_fig(gcf,plotname);
end

%% 
fprintf('--------%s----------------------------------\n',AnaFlag)
%fprintf('stat. only mnu^2 = %.3f (-%.3f + %.3f) eV^2 \n',FitResultStat.par(1),-FitResultStat.errNeg(1),FitResultStat.errPos(1));
fprintf('stat. only  \n mnu^2 = %.3f (-%.3f + %.3f) eV^2 \n E0     = %.3f eV +- %.3f \n chi2 = %.2f (%.0f dof), p-value = %.2f\n',...
    FitResultStat.par(1),-FitResultStat.errNeg(1),FitResultStat.errPos(1),...
    FitResultStat.par(2)+dCM.A.ModelObj.Q_i,FitResultStat.err(2),...
    FitResultStat.chi2min,FitResultStat.dof,1-chi2cdf(FitResultStat.chi2min,FitResultStat.dof));

fprintf('syst. only mnu^2 = %.3f (-%.3f + %.3f) eV^2 \n',FitResultCM.par(1),sqrt(FitResultCM.errNeg(1)^2-FitResultStat.errNeg(1)^2),sqrt(FitResultCM.errPos(1)^2-FitResultStat.errPos(1)^2));
fprintf('total  \n mnu^2 = %.3f (-%.3f + %.3f) eV^2 \n E0     = %.3f eV +- %.3f \n chi2 = %.2f (%.0f dof), p-value = %.2f\n',...
    FitResultCM.par(1),-FitResultCM.errNeg(1),FitResultCM.errPos(1),...
    FitResultCM.par(2)+dCM.A.ModelObj.Q_i,FitResultCM.err(2),...
    FitResultCM.chi2min,FitResultCM.dof,1-chi2cdf(FitResultCM.chi2min,FitResultCM.dof));
%%
%dCM.A.PlotFitMultiRing('PlotPar','Norm','Blind','ON','savePlot','ON');%,'Bkg'
%dCM.A.PlotFitMultiRing('PlotPar','Bkg','Blind','ON','savePlot','ON');%,'Bkg'
