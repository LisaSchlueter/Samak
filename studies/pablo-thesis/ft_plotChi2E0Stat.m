close all;
addpath(genpath('../../../Samak2.0'));

RunList = load('RunList100Good.mat');
RunList = RunList.RunList100Good;
RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];

if ~exist('SPAll','var')
    BuildTable(1);
    load('BuildTableWorkspace')
end

figchi2 = figure(8);
set(figchi2,'units','normalized','position',[0.001 0.05 0.9 0.6]);


subplot(1,3,1)
chi2_s = SPAll.TableShortStat(:,9); dof_s = SPAll.TableShortStat(1,10);
hold on
[~,Nbins,Wbins] = nhist(chi2_s);
binWidth = Wbins(2) - Wbins(1);
x_N = linspace(min(chi2_s),max(chi2_s),...
    length(Nbins)*50);
plot(x_N,length(chi2_s)*binWidth*chi2pdf(x_N,dof_s),'LineWidth',4);
hold off
xlabel('Short (E_0 - 200 eV)');
title(['DoF = ',num2str(dof_s)])

subplot(1,3,2)
chi2_s = SPAll.TableMedStat(:,9); dof_s = SPAll.TableMedStat(1,10);
hold on
[~,Nbins,Wbins] = nhist(chi2_s);
binWidth = Wbins(2) - Wbins(1);
x_N = linspace(min(chi2_s),max(chi2_s),...
    length(Nbins)*50);
plot(x_N,length(chi2_s)*binWidth*chi2pdf(x_N,dof_s),'LineWidth',4);
hold off
xlabel('Medium (E_0 - 400 eV)');
title(['DoF = ',num2str(dof_s)])


subplot(1,3,3)
chi2_s = SPAll.TableLongStat(:,9); dof_s = SPAll.TableLongStat(1,10);
hold on
[~,Nbins,Wbins] = nhist(chi2_s);
binWidth = Wbins(2) - Wbins(1);
x_N = linspace(min(chi2_s),max(chi2_s),...
    length(Nbins)*50);
plot(x_N,length(chi2_s)*binWidth*chi2pdf(x_N,dof_s),'LineWidth',4);
hold off
xlabel('Long (E_0 - 1600 eV)');
title(['DoF = ',num2str(dof_s)])

suptitle('\chi^2 Dist. Stat. All Runs');


export_fig('plots/Chi2E0Stat.pdf','-pdf')





