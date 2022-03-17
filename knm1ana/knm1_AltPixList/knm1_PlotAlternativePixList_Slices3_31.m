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
savename3 = sprintf('%sknm1_PixListAlt_%s_%s_%s_%.0feV_%s_%2g.mat',...
    savedir,'Slice3',DataType,strrep(freePar,' ',''),range,chi2,NP);
savename3_1 = sprintf('%sknm1_PixListAlt_%s_%s_%s_%.0feV_%s_%2g.mat',...
    savedir,'Slice3_1',DataType,strrep(freePar,' ',''),range,chi2,NP);



if exist(savename3,'file')
    d3 = importdata(savename3);
    fprintf('load %s \n',savename3)
    
    d31 = importdata(savename3_1);
    fprintf('load %s \n',savename3_1)
else
    fprintf('file not found %s \n',savename3)
    return
end

AngleDeg3_m = mean(d3.SliceAngPos,2);
AngleDeg31_m = mean(d31.SliceAngPos,2);
AngleDeg31_m(1) = 0;

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);

 %% m^2 as a function of pixel position
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);

xErr = max(d3.SliceAngPos')-AngleDeg3_m';

e3 = errorbar(AngleDeg3_m,d3.mNuSq,d3.mNuSqErr,d3.mNuSqErr,xErr,xErr,...
    'k.','MarkerSize',17,'CapSize',0,'LineWidth',1);
hold on
e31 = errorbar(AngleDeg31_m,d31.mNuSq,d31.mNuSqErr,d31.mNuSqErr,xErr,xErr,...
    '.','MarkerSize',17,'CapSize',0,'LineWidth',1,'Color',rgb('DimGray'));

xtickformat('%.0f°')
xlabel('Azimuthal angle');
ylabel(sprintf('{\\itm}_\\nu^2 (eV^{ 2})'));
 PrettyFigureFormat('FontSize',18);
 
 xlim([-10 360]);
 ylim([-17,10]);
 
 

% pname3 = sprintf('%sknm2_AltPixList_%s_mNuSq.pdf',pltdir,AltPixList);
% export_fig(pname3);
% fprintf('save plot to %s \n',pname3);

 