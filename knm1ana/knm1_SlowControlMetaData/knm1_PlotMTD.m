% plot MTD for knm1
% for PhD thesis
savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);
FontSize = 22;
if exist(savefile,'file')
    load(savefile);
else
    return
end

Mode = 'Rel'; % y-axis. relative time fraction or actual cumulative hours

close all
switch Mode
    case 'Abs'
        TimeHour = R.RunData.TimeSec./(60*60);
    case 'Rel'
        TimeHour = 1;
end
f1 = figure('Units','normalized','Position',[1.1,1.1,0.5,0.4]);%  GetFigure;
%b200 = bar(R.RunData.qU_RM,1e2*R.RunData.qUfrac_RM,'EdgeColor',rgb('Red'),'FaceColor',rgb('Red'),'FaceAlpha',0.8);
pE0 = plot(18574.*ones(2,1),TimeHour*linspace(0,1.5*max(R.RunData.qUfrac),2),':','LineWidth',3,'Color',rgb('Orange'));
hold on;
b90 = bar(R.RunData.qU(1:end),TimeHour*R.RunData.qUfrac(1:end),'EdgeColor',rgb('Silver'),'FaceColor',rgb('Silver'),'FaceAlpha',0.8,'BarWidth',1);
b40 = bar(R.RunData.qU(R.exclDataStart:end),TimeHour*R.RunData.qUfrac(R.exclDataStart:end),'EdgeColor',rgb('DodgerBlue'),'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.8,'BarWidth',b90.BarWidth);
bbkg = bar(R.RunData.qU(end-4:end),TimeHour*R.RunData.qUfrac(end-4:end),'EdgeColor',rgb('SkyBlue'),'FaceColor',rgb('SkyBlue'),'FaceAlpha',0.8,'BarWidth',0.196);
PrettyFigureFormat('FontSize',FontSize);
pNone90 = plot(NaN,NaN,'w.');pNone40 = plot(NaN,NaN,'w.');pNonebkg = plot(NaN,NaN,'w.');pNone = plot(NaN,NaN,'w.');
% legend
leg = legend([b90,b40,bbkg,pE0],...
    sprintf('Full interval \n[{\\itE}_0 - %.0f eV, {\\itE}_0 + %.0f eV]',abs(R.ModelObj.qU(1)-18574),R.ModelObj.qU(end)-18574),...
    sprintf('Analysis interval \n[{\\itE}_0 - %.0f eV, {\\itE}_0 + %.0f eV]',abs(R.ModelObj.qU(R.exclDataStart)-18574),R.ModelObj.qU(end)-18574),...
    sprintf('Background \n[{\\itE}_0 + %.0f eV,  {\\itE}_0 + %.0f eV]',R.ModelObj.qU(end-4)-18574,R.ModelObj.qU(end)-18574),...
    sprintf('{\\itE}_0 = %.0f eV',18574),...
      'Location','west');
legend boxoff
leg.NumColumns = 1;
% style
xlabel('Retarding energy (eV)');
switch Mode
    case 'Abs'
      ylabel('Measurement time (hours)') 
      ylim([0 45]);
    case 'Rel'
       ylabel('Relative measurement time')
       ylim([0 0.08]);
end


xlim([18477 18625])
ax = gca;
ax.XAxis.Exponent = 0;

leg.Position(4) = 0.45;
leg.Position(2) = 0.33;
leg.ItemTokenSize(2) = 10;
set(gca,'XMinorTick','off');


%%
plotdir =  [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/plots/'];
MakeDir(plotdir)

switch Mode
    case 'Abs'
      plotfile = sprintf('%sknm1_MTD.pdf',plotdir);
    case 'Rel'
     plotfile = sprintf('%sknm1_MTD_rel.pdf',plotdir);
end

export_fig(plotfile);
fprintf('save plot to %s \n',plotfile);

