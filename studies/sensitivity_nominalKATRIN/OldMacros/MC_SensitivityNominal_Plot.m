%--------------------------------------------------------------------------------------
%  Sensitivity Monte Carlo Study 
% Nominal KATRIN settings
% Plot Fit Distributions and calculate sensitivity
%--------------------------------------------------------------------------------------
MACE_Ba_T = 7*1e-04;%1e-04.*(3:10);
WGTS_B_T = 0.7*3.6;
%Sensitivity90 = zeros(numel(MACE_Ba_T_all),1);
%for i=1:numel(MACE_Ba_T_all)
%MACE_Ba_T = MACE_Ba_T_all(i);   % B-Field analyzing plane. Background is scaled automatically in SensitivityStudy_NominalKATRIN
range = 30;                 % MTD energy range: 30,45,60eV below E0 (18575eV)
chi2 = 'chi2Stat';
SysEffect ='FSD';
SysFluct = 'OFF';
if strcmp(chi2,'chi2Stat')
    SysEffect = '';
    SysFluct = 'OFF';
end
nSamples = 10000;%5000;
TimeSec = 3*365*24*60*60;
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
file_name = sprintf('./../sensitivity_nominalKATRIN/results/SensitivityStudy_NominalKATRIN_MTD-%s_Time-%.2fy_mnu%.1feV_StatFluct%s_SysFluct%s_%.0fSamples_%s%s.mat',...
        TD,TimeSec/(60*60*24*365),0,'ON',SysFluct,nSamples,chi2,SysEffect);
try
    load(file_name,'-mat');
catch
    fprintf('File doesnt exist!\n')
    return
end
BKG_RateSim = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T);
mNu = sort(par(1,:)); % Neutrino mass^2
mNu90 = mNu(nSamples*0.9); %neutrino mass squared sensitivity
mNu90low = mNu(nSamples*0.05);
mNu90up  = mNu(nSamples*0.95);
Sensitivity90(1) = sqrt(mNu90);

if strcmp(chi2,'chi2Stat')
    SysEffect = '';
    chi2label = 'stat';
else
    chi2label = 'stat + sys';
end
switch SysFluct
    case 'ON'
        maintitle =  sprintf('Samak %s Fit to Simulation (asimov + stat. + sys. fluct) \n - KATRIN %.0f years - Scan Range %.0f eV - Ba %.0fG - BKG %.0fmcps - %.0f Samples -',...
            chi2label,TimeSec/(60*60*24*365),range,MACE_Ba_T*1e4,BKG_RateSim*1e3,nSamples);
    case 'OFF'
        maintitle =  sprintf('Samak %s Fit to Simulation (asimov + stat. fluct) \n - KATRIN %.0f years - Scan Range %.0f eV - Ba %.0fG  - BKG %.0fmcps - %.0f Samples -',...
            chi2label,TimeSec/(60*60*24*365),range,MACE_Ba_T*1e4,BKG_RateSim*1e3,nSamples);
end
%%
Mode = 'TwoSided';%'OneSided';
Legend = 'OFF';
BinWidth = 0.008;
PlotFontSize = 18;
f10 = figure('Name','NuMass','Renderer','opengl');
%f10 = figure(10)
set(f10, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
h1 = histogram(mNu,'FaceColor',rgb('CadetBlue'),'BinWidth',BinWidth);
Counts = h1(1).get.BinCounts;
hold on;
switch Mode
    case 'OneSided'
        % find bin index, where cumulative sum = 90% (for Display)
        Index90 = find(abs(h1.BinEdges-mNu90) == min(abs(h1.BinEdges-mNu90)));
        mNu90hist = h1.BinEdges(Index90); %neutrino mass^2 90% from histogram
        h2 = histogram(mNu(mNu<=mNu90hist),'FaceColor',rgb('DarkRed'),'BinWidth',BinWidth);
        %x1 = [0.675, 0.675];
        x1 = ((mNu90hist+0.3)/0.6)*[1 1]-0.5*BinWidth-0.004;
        y1 = [0.8, Counts(Index90-1)/max(ylim)+0.02];
        a1= annotation(f10,'textarrow',x1,y1,'String',sprintf('c'));
        a1.FontSize = PlotFontSize+6;
        a1.FontWeight = 'bold';
        a1.LineWidth = 1.5;
    case 'TwoSided' 
        Index90up = find(abs(h1.BinEdges-mNu90up) == min(abs(h1.BinEdges-mNu90up)));
        Index90low = find(abs(h1.BinEdges-mNu90low) == min(abs(h1.BinEdges-mNu90low)));
        mNu90histup = h1.BinEdges(Index90up); %neutrino mass^2 90% from histogram
        mNu90histlow = h1.BinEdges(Index90low); %neutrino mass^2 90% from histogram
        mNuUp_tmp = mNu(mNu<=mNu90histup);
        h2 = histogram(mNuUp_tmp(mNuUp_tmp>=mNu90histlow),'FaceColor',rgb('DarkRed'),'BinWidth',BinWidth);
        xLow = [0.395, 0.395];
        %xLow = (abs(mNu90histlow)+0.3)/0.6*[1 1]+0.13-0.5; % upper
        yLow = [0.8, Counts(Index90low-2)/max(ylim)+0.15];
        xUp = ((mNu90histup+0.3)/0.6)*[1 1]-2*BinWidth; % upper
        yUp = [0.8, Counts(Index90up)/max(ylim)+0.08];
        aLow= annotation(f10,'textarrow',xLow,yLow,'String',sprintf('a'),'Units','normalized');
        aUp= annotation(f10,'textarrow',xUp,yUp,'String',sprintf('b'),'Units','normalized');
        aUp.FontSize = PlotFontSize+6; aLow.FontSize = PlotFontSize+6;
        aUp.FontWeight = 'bold';  aLow.FontWeight = 'bold';
        aUp.LineWidth = 1.5;  aLow.LineWidth = 1.5;
end
%plot([0 0],[0 45],'-','LineWidth',2,'Color',rgb('DarkSlateGray'));
%plot([mNu90low mNu90low]-0.5*BinWidth,[0 40],'--','LineWidth',2);
%plot([mNu90up mNu90up]-0.5*BinWidth,[0 40],'--','LineWidth',2);
x= [min(mNu)-0.05:0.01:max(mNu)+0.05];
pnorm = plot(x+0.5*BinWidth,nSamples*BinWidth.*normpdf(x,mean(mNu),std(mNu)),'Color',rgb('Black'),'LineWidth',3);
PrettyFigureFormat;
set(gca,'FontSize',PlotFontSize);
xlabel('m_{\nu}^2 (eV^2)');
ylabel('number of samples');
switch Legend
    case 'ON'
        leg1 = ['Distribution <m_{\nu}^2>=',sprintf('%.3f eV^2 \\pm %.3f eV^2 (std)',mean(mNu),std(mNu))];
        leg2 = ['\sigma(m_{\nu}^2)=',sprintf('%.3f eV² (90%% C.L.) \n',mNu90),'\sigma(m_{\nu})=',sprintf('%.3f eV (90%% C.L.) \n',sqrt(mNu90))] ;
        legend([h1,h2],leg1,leg2)
        title(maintitle);
    case 'OFF'
        leg = legend(h2,'90% C.L.','Location','best');
        legend boxoff;
end

save_name = sprintf('SensitivityMC_%s_%s%s_%s_%.0f-BKG_%.2fyears_%.0fSamples',Mode,chi2,SysEffect, TD,BKG_RateSim*1e3,TimeSec/(60*60*24*365),nSamples);
export_fig(f10,['./plots/png/MC/',save_name,'.png']);
savefig(f10,['./plots/fig/MC/',save_name,'.fig'],'compact');
publish_figurePDF(f10,['./plots/pdf/MC/',save_name,'.pdf']);
%% Test of distribution is gaussian
% chi2 Goodness of Fit
[h,p, stats] = chi2gof(sort(mNu),'NBins',nSamples/10);
if h==0
    fprintf('nu-mass² distribution is gaussian with significane of %.2g %% (p-value) \n',p*100)
else
   fprintf('nu-mass² distribution is NOT gaussian with significane of %.2g %% (1 - p-value) \n',(1-p)*100)
end
%% Plot other fit paramete
f11 = figure(11);
set(f11, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';

subplot(2,4,[1,2]);
h1 = histfit(mNu,nSamples/100,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
leg = legend(sprintf('<m_{\\nu}^2 > = %.2g \\pm %.1g (std) eV^2',mean(par(1,:)),std(par(1,:))),'normal');
leg.Location = 'best'; leg.FontSize = PlotFontSize-3;
legend boxoff
xlabel('m_{\nu}^2 (eV^2)');
PrettyFigureFormat;
set(gca,'FontSize',PlotFontSize);

subplot(2,4,[3,4]);
h1 = histfit(par(2,:)+18573.7,nSamples/100,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
xlabel('E_{0eff} (eV)');
leg = legend(sprintf('<E_{0eff}> = %.1f \\pm %.1g (std) eV',mean(par(2,:)+18573.7),std(par(2,:))),'normal');
leg.Location = 'best'; leg.FontSize = PlotFontSize-3;legend boxoff;
PrettyFigureFormat;
set(gca,'FontSize',PlotFontSize);

subplot(2,4,5);
h1 = histfit(par(3,:)+BKG_RateSim,nSamples/100,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
xlabel('BKG (cps)');
leg = legend(sprintf('<BKG> = %.1g \\pm %.1g (std) cps',mean(par(3,:)+BKG_RateSim),std(par(3,:))),'normal');
leg.Location = 'best'; leg.FontSize = PlotFontSize-3;legend boxoff;
PrettyFigureFormat;
set(gca,'FontSize',PlotFontSize);

subplot(2,4,6);
h1 = histfit(par(4,:)+1,nSamples/100,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
xlabel('N');
leg = legend(sprintf('<N> = %.1f \\pm %.1g (std)',mean(par(4,:)+1),std(par(4,:))),'normal');
leg.Location = 'best'; leg.FontSize = PlotFontSize-3;legend boxoff;
PrettyFigureFormat;
set(gca,'FontSize',PlotFontSize);

subplot(2,4,[7,8])
h1 = histogram(chi2min,'FaceColor',rgb('CadetBlue'));
hold on
xlabel(sprintf('\\chi^2 (%.0f dof)',dof(1)));
PrettyFigureFormat;
x =(min(chi2min):max(chi2min))';
y = repmat(dof(1),numel(x),1);
pchi2 = plot(x,chi2pdf(x,y)*h1.BinWidth*nSamples,'LineWidth',3,'Color',rgb('Goldenrod'));
leg = legend([h1, pchi2],...
    sprintf('<\\chi2> = %.1f \\pm %.1g (std)',mean(chi2min),std(chi2min)),...
    sprintf('\\chi2 pdf for %.0f dof',dof(1)));
leg.Location = 'best'; leg.FontSize = PlotFontSize-3;
legend boxoff;
set(gca,'FontSize',PlotFontSize);

save_name = sprintf('SensitivityMC_ParOverview_%s%s_%s_%.0f-BKG_%.2fyears_%.0fSamples',chi2,SysEffect, TD,BKG_RateSim*1e3,TimeSec/(60*60*24*365),nSamples);
export_fig(f11,['./plots/png/',save_name,'.png']);
savefig(f11,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f11,['./plots/pdf/',save_name,'.pdf']);
%end
%%
% figure(123);
% BKG_all = GetBackground(MACE_Ba_T);
% x= BKG_all*1e3;% MACE_Ba_T_all;
% plot(x,Sensitivity90*1e3,'--x','LineWidth',3,'MarkerSize',8);
% %xlabel('B_a (10^{-4} T)');
% xlabel('Background (mcps)');
% ylabel('\sigma(m_{\nu}) 90% C.L. (meV)')
% PrettyFigureFormat;
% set(gca,'FontSize',16);
% xlim([min(x), max(x)]);
% grid on;