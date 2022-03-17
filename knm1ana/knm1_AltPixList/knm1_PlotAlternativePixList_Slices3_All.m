% plot spin-off from to knm1_PlotAlternativePixList_Slices 
% compare slices 3 and slice 3_1

freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;                % fit range in eV below endpoint
chi2 = 'chi2Stat';
NP = 1.064;
RecomputeFlag = 'OFF';

% label
savedir = [getenv('SamakPath'),'knm1ana/knm1_AltPixList/results/'];
savename = sprintf('%sknm1_PixListAlt_%s_%s_%s_%.0feV_%s_%2g.mat',...
    savedir,'Slice',DataType,strrep(freePar,' ',''),range,chi2,NP);
savename3_1 = sprintf('%sknm1_PixListAlt_%s_%s_%s_%.0feV_%s_%2g.mat',...
    savedir,'Slice3_1',DataType,strrep(freePar,' ',''),range,chi2,NP);



if exist(savename,'file')
    d = importdata(savename);
    fprintf('load %s \n',savename)
    
    d31 = importdata(savename3_1);
    fprintf('load %s \n',savename3_1)
else
    fprintf('file not found %s \n',savename)
    return
end

AngleDeg_m = mean(d.SliceAngPos,2);
AngleDeg31_m = mean(d31.SliceAngPos,2);
AngleDeg31_m(1) = 0;

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);

if any(d.mNuSqErr==0)
    %sometimes asymmetric uncertainty may not work properly (convergence problems...), use some average uncertainty....
    fprintf(2,'WARNING %.0f slices have convergence problems (->tiny uncertainty) \n',sum(d.mNuSqErr==0));
    d.mNuSqErr(d.mNuSqErr==0)= median(d.FitResult.err(:,1));
end
% exlucde slices that didn't converge
InclIdx = d.mNuSqErr<2.5*median(d.mNuSqErr); %& d.mNuSqErr>0;
%% m^2 as a function of pixel position
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);


eall = errorbar(AngleDeg_m(InclIdx),d.mNuSq(InclIdx),d.mNuSqErr(InclIdx),...
    'k.','MarkerSize',17,'CapSize',0,'LineWidth',1,'Color',rgb('Silver'));
hold on
e31 = errorbar(AngleDeg31_m,d31.mNuSq,d31.mNuSqErr,...
    '.','MarkerSize',17,'CapSize',0,'LineWidth',1,'Color',rgb('DodgerBlue'));

xtickformat('%.0f°')
xlabel('Azimuthal angle');
ylabel(sprintf('{\\itm}_\\nu^2 (eV^{ 2})'));
PrettyFigureFormat('FontSize',18);

xticks([0:45:360]);
xlim([-10 360]);
%ylim([-17,10]);

if any(InclIdx==0)
    angle = AngleDeg_m(~InclIdx);
    mnu = d.mNuSq(~InclIdx);
    mnuerr = d.mNuSqErr(~InclIdx);
    ym = ylim;
    for i=1:numel(mnu)
          t = text(5,ym(1)+5*i,sprintf('Slice at \\langle%.0f^\\circ\\rangle: {\\itm}^2_\\nu = (%.0f\\pm%.0f) eV^2',angle(i),mnu(i),mnuerr(i)),...
            'FontSize',16,'Color',rgb('Silver'));
    end
end

% pname3 = sprintf('%sknm2_AltPixList_%s_mNuSq.pdf',pltdir,AltPixList);
% export_fig(pname3);
% fprintf('save plot to %s \n',pname3);

 