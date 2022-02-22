% plot KNM2 MTD
Mode = 'Rel';
FontSize = 22;

% get Data
savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s_SysBudget40.mat',...
    savedir,3,'Real',40,'mNuE0BkgNorm','chi2CMShape','StackPixel','KNM2_0p1eV');
d = importdata(savename);

% plot
close all
switch Mode
    case 'Abs'
        TimeHour = d.A.RunData.TimeSec./(60*60);
    case 'Rel'
        TimeHour = 1;
end
f1 = figure('Units','normalized','Position',[1.1,1.1,0.5,0.4]);%  GetFigure;
%b200 = bar(d.A.RunData.qU_RM,1e2*d.A.RunData.qUfrac_RM,'EdgeColor',rgb('Red'),'FaceColor',rgb('Red'),'FaceAlpha',0.8);
pE0 = plot(18574.*ones(2,1),TimeHour*linspace(0,1.5*max(d.A.RunData.qUfrac),2),':','LineWidth',3,'Color',rgb('Orange'));
hold on;
b90 = bar(d.A.RunData.qU(1:end),TimeHour*d.A.RunData.qUfrac(1:end),'EdgeColor',rgb('Silver'),'FaceColor',rgb('Silver'),'FaceAlpha',0.8,'BarWidth',1);
b40 = bar(d.A.RunData.qU(d.A.exclDataStart:end),TimeHour*d.A.RunData.qUfrac(d.A.exclDataStart:end),'EdgeColor',rgb('DodgerBlue'),'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.8,'BarWidth',b90.BarWidth);
bbkg = bar(d.A.RunData.qU(end-4:end),TimeHour*d.A.RunData.qUfrac(end-4:end),'EdgeColor',rgb('SkyBlue'),'FaceColor',rgb('SkyBlue'),'FaceAlpha',0.8,'BarWidth',0.196);
PrettyFigureFormat('FontSize',FontSize);
pNone90 = plot(NaN,NaN,'w.');pNone40 = plot(NaN,NaN,'w.');pNonebkg = plot(NaN,NaN,'w.');pNone = plot(NaN,NaN,'w.');
% legend
leg = legend([b90,b40,bbkg,pE0],...
    sprintf('Full interval \n[{\\itE}_0 - %.0f eV, {\\itE}_0 + %.0f eV]',abs(d.A.ModelObj.qU(1)-18574),d.A.ModelObj.qU(end)-18574),...
    sprintf('Analysis interval \n[{\\itE}_0 - %.0f eV, {\\itE}_0 + %.0f eV]',abs(d.A.ModelObj.qU(d.A.exclDataStart)-18574),d.A.ModelObj.qU(end)-18574),...
    sprintf('Background \n[{\\itE}_0 + %.0f eV,  {\\itE}_0 + %.0f eV]',d.A.ModelObj.qU(end-4)-18574,d.A.ModelObj.qU(end)-18574),...
    sprintf('{\\itE}_0 = %.0f eV',18574),...
      'Location','west');
legend boxoff
leg.NumColumns = 1;
% style
xlabel('Retarding energy (eV)');
switch Mode
    case 'Abs'
      ylabel('Measurement time (hours)') 
      ylim([0 67]);
    case 'Rel'
       ylabel('Relative measurement time')
       ylim([0 0.08]);
end


xlim([18447 18574+140])
ax = gca;
ax.XAxis.Exponent = 0;
xticks([18480:50:19000])
leg.Position(4) = 0.45;
leg.Position(2) = 0.33;
leg.ItemTokenSize(2) = 10;
set(gca,'XMinorTick','off');


%%
plotdir =  [getenv('SamakPath'),'knm2ana/knm2_DataSet/plots/'];
MakeDir(plotdir)

switch Mode
    case 'Abs'
      plotfile = sprintf('%sknm2_MTD.pdf',plotdir);
    case 'Rel'
     plotfile = sprintf('%sknm2_MTD_rel.pdf',plotdir);
end

export_fig(plotfile);
fprintf('save plot to %s \n',plotfile);


%% info:
fprintf('%.0f scan-steps (full), %.0f analysis range \n',numel(d.A.RunData.qU),numel(d.A.RunData.qU(d.A.exclDataStart:end)));
fprintf('time above endpoint %.1f%% (full), %.1f%% (analysis interval) \n',...
    1e2.*sum(d.A.RunData.qUfrac(d.A.RunData.qU>18574))./sum(d.A.RunData.qUfrac),...
     1e2.*sum(d.A.RunData.qUfrac(d.A.RunData.qU>18574))./sum(d.A.RunData.qUfrac(d.A.exclDataStart:end)));
