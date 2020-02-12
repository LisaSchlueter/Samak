%%Samak
range = 90;
Rings = 1:12;
Sfilepath = sprintf('../KNM1_RingAnalysis/results/');
Sfilename = sprintf('FitRings175mV_%.0frange.mat',range);
s = importdata([Sfilepath,Sfilename]);
E0Samak    = s.E0-mean(s.E0);
E0ErrSamak =s.E0Err;

[Spar, Serr, Schi2min,Sdof] =linFit(Rings',E0Samak,E0ErrSamak);
%% Fitritum
Ffilename = sprintf('fitrium_ringwise_%.0f.txt',range);
Ffilepath = [Sfilepath,'/FitriumRingwiseResults/'];
f = importdata([Ffilepath,Ffilename]);
E0Fitrium = f.data(:,2);
E0ErrFitrium = f.data(:,3);

[Fpar, Ferr, Fchi2min,Fdof] =linFit(Rings',E0Fitrium,E0ErrFitrium);
%% KaFit
Kfilename = sprintf('KaFitResult_Run51279_51326_qU%.0feV.txt',range);
Kfilepath = [Sfilepath,'/KaFitResult/'];
k = importdata([Kfilepath,Kfilename]);
E0kafit = k.data(:,2)-mean(k.data(:,2));
E0Errkafit = k.data(:,3);
[Kpar, Kerr, Kchi2min,Kdof] =linFit(Rings',E0kafit,E0Errkafit);
%%

fig1 = figure('Renderer','opengl');
set(fig1,'units','normalized','pos',[0.1, 0.1,0.7,0.7]);
CommonArg = {'o','LineWidth',3,'LineStyle','none','MarkerSize',7};
pk = errorbar(Rings,E0kafit,E0Errkafit,CommonArg{:},...
   'Color',rgb('SkyBlue'),'MarkerFaceColor',rgb('SkyBlue'));
hold on;
lk = plot(Rings,Kpar(1).*Rings+Kpar(2),'--','Color',rgb('SkyBlue'),'LineWidth',2);
pf = errorbar(Rings-0.15,E0Fitrium,E0ErrFitrium,CommonArg{:},...
     'Color',rgb('Orange'),'MarkerFaceColor',rgb('Orange'));
lf = plot(Rings,Fpar(1).*Rings+Fpar(2),'--','Color',rgb('Orange'),'LineWidth',2);
ps = errorbar(Rings+0.15,E0Samak,E0ErrSamak,CommonArg{:},...
     'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('SlateGray'));
ls = plot(Rings,Spar(1).*Rings+Spar(2),'--','Color',rgb('SlateGray'),'LineWidth',2);
hold off;
leg = legend([pk,pf,ps,lk,lf,ls],'KaFit','Fitrium','Samak',...
    sprintf('slope: (%.1f \\pm %.1f) meV  %.1f /%.0f dof',Kpar(1)*1e3,Kerr(1)*1e3,Kchi2min,Kdof),...
    sprintf('slope: (%.1f \\pm %.1f) meV  %.1f /%.0f dof',Fpar(1)*1e3,Ferr(1)*1e3,Fchi2min,Fdof),...
    sprintf('slope: (%.1f \\pm %.1f) meV  %.1f /%.0f dof',Spar(1)*1e3,Serr(1)*1e3,Schi2min,Sdof));
legend boxoff
leg.Location = 'northwest';
leg.NumColumns = 2;
xlim([0.5,12.5])
ylabel('E_0 - <E_0> (eV)')
xlabel('ring')
PrettyFigureFormat;
title(sprintf('%.0f runs (%.0f - %.0f) - %.0f pixels total - %.0feV range',numel(s.M(1).StackedRuns),s.M(1).RunList(1),s.M(1).RunList(end),118,range));
set(gca,'FontSize',20);

print(sprintf('../KNM1_RingAnalysis/plots/FitterComparison_RingFit_175mVRuns_%.0feVrange.png',range),'-dpng','-r500');
