%--------------------------------------------------------------------------------------
%  Compute and Plot Sensitivity from Monte Carlo Method
%  Plot Parameter Distribution
%  First Tritium KATRIN settings
%--------------------------------------------------------------------------------------
%SysEffect_all  = {'TC','TASR','FSD','RF','all'};
%for i=1:numel(SysEffect_all)
 %  SysEffect = SysEffect_all{i};
 SysEffect = 'all';
nSamples = 5000;
RunList = 'StackCD100all';
StatFluct = 'ON';
SysFluct = 'ON';
chi2 = 'chi2Stat';
SysEffect  = '';%'TC','TASR','FSD','RF','all'};
fixPar = '1 5 6';
Q_i    = 18573.7;
exclDataStart = 7;

switch exclDataStart
    case 7
        belowE0 = 402;
    case 9
        belowE0 = 202;
end
switch chi2
    case 'chi2Stat'
        SysEffect = '';
        SysFluct = 'OFF';
        chi2label = 'stat';
    case 'chi2CM'
        chi2label = sprintf('stat + sys (%s CM)',SysEffect);
end
%---------------------------------inputs end-----------------------------------
%Load Results
file_name = sprintf('./results/MC_SensitivityStudy_FTKATRIN_%s_StatFluct%s_SysFluct%s_%s%s_fixPar%s_%.0feVrange_Qi%.0f_%.0fSamples.mat',...
     RunList,StatFluct,SysFluct,chi2,SysEffect,strrep(fixPar,' ',''),belowE0,Q_i*10,nSamples);
 
try
    load(file_name,'-mat');
    E0      = sort(cellfun(@(x) x.par(2),FitResult)+Q_i);
    E0Err   = cellfun(@(x) x.err(2),FitResult);
    N       = cellfun(@(x) x.par(4)+1,FitResult);
    BKG     = cellfun(@(x) x.par(3),FitResult);
    dof     = cellfun(@(x) x.dof,FitResult);
    chi2min = cellfun(@(x) x.chi2min,FitResult);
catch
    fprintf('File doesnt exist!\n')
    return
end

% Compute Sensitivity
E090low = abs(E0(nSamples*0.05)-mean(E0)); % Sensitivity 90% C.L.
E090up  = E0(nSamples*0.95)-mean(E0);
%% Plot E0 distribution
switch SysFluct
    case 'ON'
        maintitle =  sprintf('Samak %s Fit to Simulation (asimov + stat. + sys. fluct) \n - KATRIN First Tritium - MTD: %s - Scan Range %.0f eV - %.0f Samples -',...
            chi2label,RunList,belowE0,nSamples);
    case 'OFF'
       maintitle =  sprintf('Samak %s Fit to Simulation (asimov + stat. fluct) \n - KATRIN First Tritium - MTD: %s - Scan Range %.0f eV - %.0f Samples -',...
            chi2label,RunList,belowE0,nSamples);
end
close all;
f10 = figure(10);
set(f10, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
h1 = histogram(E0,'FaceColor',rgb('CadetBlue'),'BinWidth',0.02);
% find bin index, where cumulative sum = 5% and 95% (for Display)
Index90low = find(abs(0.05- cumsum(h1.Values)./sum(h1.Values))==min(abs(0.05- cumsum(h1.Values)./sum(h1.Values))));
Index90up = find(abs(0.95- cumsum(h1.Values)./sum(h1.Values))==min(abs(0.95- cumsum(h1.Values)./sum(h1.Values))));
E090histlow = h1.BinEdges(Index90low); %neutrino mass^2 90%
E090histup = h1.BinEdges(Index90up);
hold on;
E0up = E0(E0<=E090histup);
h2 = histogram(E0up(E0up>=E090histlow),'FaceColor',rgb('IndianRed'),'BinWidth',0.02);
leg1 = ['Distribution <E_{0eff}> =',sprintf('%.3f eV \\pm %.3f eV (std)',mean(E0),std(E0))];
leg2 = ['\sigma(E_{0eff})=',sprintf('- %.3f + %.3f eV (90%% C.L.) \n',E090low, E090up)] ;
legend([h1,h2],leg1,leg2)
PrettyFigureFormat;
set(gca,'FontSize',16);
xlabel('E_{0eff} (eV)');
title(maintitle);

save_name = sprintf('MC_SensitivityFT_Endpoint_%s_%s%s_%.0feVrange_%.0fSamples',RunList,chi2,SysEffect,belowE0,nSamples);
export_fig(f10,['./plots/png/',save_name,'.png']);
savefig(f10,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f10,['./plots/pdf/',save_name,'.pdf']);
%% Plot Distrbution of all Fit Parameters
f11 = figure(11);
set(f11, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=18;a.FontWeight='bold';
binfactor = 50;
subplot(2,4,[1,2]); % Endpoint
h1 = histfit(E0,nSamples/binfactor,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
legend(sprintf('<E_{0eff}> = %.1f \\pm %.1g (std) eV',mean(E0),std(E0)),'normal');
xlabel('E_{0eff} (eV)');
PrettyFigureFormat;

subplot(2,4,[3,4]); %Background Bias
h1 = histfit(BKG*1e3,nSamples/binfactor,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
xlabel('Background Bias (mcps)');
legend(sprintf('<B> = %.1f \\pm %.1g (std) eV',mean(BKG),std(BKG)),'normal');
PrettyFigureFormat;

subplot(2,4,[5,6]); %Normalization 
h1 = histfit(N,nSamples/binfactor,'normal');
h1(1).FaceColor = rgb('CadetBlue');
h1(2).Color = rgb('Goldenrod'); h1(2).LineWidth = 3;
xlabel('Normalization');
legend(sprintf('<N> = %.1g \\pm %.1g (std) cps',mean(N),std(N)),'normal');
PrettyFigureFormat;

subplot(2,4,[7,8]) %chi2
h1 = histogram(chi2min,'FaceColor',rgb('CadetBlue'));
hold on
xlabel(sprintf('\\chi^2 (%.0f dof)',dof(1)));
PrettyFigureFormat;
x =(min(chi2min):max(chi2min))';
y = repmat(dof(1),numel(x),1);
pchi2 = plot(x,chi2pdf(x,y)*h1.BinWidth*nSamples,'LineWidth',3,'Color',rgb('Goldenrod'));
legend([h1, pchi2],...
    sprintf('<\\chi2> = %.1f \\pm %.1g (std)',mean(chi2min),std(chi2min)),...
    sprintf('\\chi2 pdf for %.0f dof',dof(1)));

save_name = sprintf('MC_SensitivityFT_ParOverview_%s_%s%s_%.0feVrange_%.0fSamples',RunList,chi2,SysEffect,belowE0,nSamples);
export_fig(f11,['./plots/png/',save_name,'.png']);
savefig(f11,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f11,['./plots/pdf/',save_name,'.pdf']);

%end



