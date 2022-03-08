% p

freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;                % fit range in eV below endpoint
chi2 = 'chi2Stat';
NP = 1.064;

savedir = [getenv('SamakPath'),'knm1ana/knm1_AltPixList/results/'];
savename = sprintf('%sknm1_SinglePixFit_%s_%s_%.0feV_%s_%2g.mat',...
    savedir,DataType,strrep(freePar,' ',''),range,chi2,NP);

if exist(savename,'file') 
    d = importdata(savename);
    fprintf('load %s \n',savename)
else
    fprintf('file not found %s \n',savename)
    return
end
%%
%plot(d.mNuSq)
FitResults = d.FitResult(1:numel(d.PixList));
mNuSqErr = cell2mat(cellfun(@(x) x.err(1),FitResults,'UniformOutput',false));
mNuSq = cell2mat(cellfun(@(x) x.par(1),FitResults,'UniformOutput',false));
histogram(mNuSqErr)


%%
close all
x = d.mNuSq;
x(abs(x)>70) = NaN;
[plotHandle, cbHandle,AngleRad] = FPDViewer(x,'Label','ON','ReDrawSkeleton','ON');
%% exclude some  pixels, convergence issues
 InclIdx = find(mNuSqErr>3 & mNuSqErr<50);
 
 mNuSq = mNuSq(InclIdx);
 mNuSqErr = mNuSqErr(InclIdx);
   PixList = d.PixList(InclIdx);
 %%
 AngleDeg = rad2deg(AngleRad);
 AngleDeg = AngleDeg+(45-AngleDeg(1));
 close all
 GetFigure;
 AngleDeg_plt = AngleDeg(PixList);
 
 
 
 a1 = area([0 90],[100 100],'BaseValue',-200,'FaceColor',rgb('Orange'),'FaceAlpha',0.5,'EdgeColor','none');
 hold on;
 a2 = area([90 180],[100 100],'BaseValue',-200,'FaceColor',rgb('FireBrick'),'FaceAlpha',0.5,'EdgeColor','none');
 a3 = area([180 270],[100 100],'BaseValue',-200,'FaceColor',rgb('DodgerBlue'),'FaceAlpha',0.5,'EdgeColor','none');
 a4 = area([270 360],[100 100],'BaseValue',-200,'FaceColor',rgb('ForestGreen'),'FaceAlpha',0.5,'EdgeColor','none');
 
 wmeanAll = wmean(mNuSq,1./mNuSqErr.^2);
 wmean1 = wmean(mNuSq(AngleDeg_plt<=90),1./mNuSqErr(AngleDeg_plt<=90).^2);
 wmean2 =wmean(mNuSq(AngleDeg_plt>90 & AngleDeg_plt<=180),1./mNuSqErr(AngleDeg_plt>90 & AngleDeg_plt<=180).^2);
 wmean3 = wmean(mNuSq(AngleDeg_plt>180 & AngleDeg_plt<=270),1./mNuSqErr(AngleDeg_plt>180 & AngleDeg_plt<=270).^2);
 wmean4 = wmean(mNuSq(AngleDeg_plt>270 & AngleDeg_plt<=360),1./mNuSqErr(AngleDeg_plt>270 & AngleDeg_plt<=360).^2);
 
 plot([0 90],wmean1.*ones(2,1),'LineWidth',3,'Color',a1.FaceColor);
 plot([90 180],wmean2.*ones(2,1),'LineWidth',3,'Color',a2.FaceColor);
 plot([180 270],wmean3.*ones(2,1),'LineWidth',3,'Color',a3.FaceColor);
 plot([270 360],wmean4.*ones(2,1),'LineWidth',3,'Color',a4.FaceColor);
 
 e1 = errorbar(AngleDeg_plt,mNuSq,mNuSqErr,'k.','MarkerSize',1,'CapSize',0,'LineWidth',1,'MarkerEdgeColor','none','MarkerFaceColor','none');
 scatter(AngleDeg(PixList),mNuSq,'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
 xlabel('Angular pixel position (degree)');
 ylabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
 PrettyFigureFormat;
 text(45,62,sprintf('Northeast \n\\langle{\\itm}_\\nu^2\\rangle = %.1f eV^2',wmean1),'HorizontalAlignment','center','FontSize',get(gca,'FontSize'));
 text(135,62,sprintf('Southeast \n\\langle{\\itm}_\\nu^2\\rangle = %.1f eV^2',wmean2),'HorizontalAlignment','center','FontSize',get(gca,'FontSize'));
 text(225,62,sprintf('Southwest \n\\langle{\\itm}_\\nu^2\\rangle = %.1f eV^2',wmean3),'HorizontalAlignment','center','FontSize',get(gca,'FontSize'));
 text(315,62,sprintf('Northwest \n\\langle{\\itm}_\\nu^2\\rangle = %.1f eV^2',wmean4),'HorizontalAlignment','center','FontSize',get(gca,'FontSize'));
 
xlim([-10 370]);
ylim([-140,82]);
