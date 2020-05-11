
%% settings
CL = 95;
range =  65;%
nGridSteps = 50;
chi2Str = 'chi2Stat';%CMShape';
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
pullFlag = 99;

[mnu4Sq,sin2T4,chi2,chi2_ref,savefileGrid] = KSN1GridSearch('range',range,...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str,...
    'DataType',DataType,...
    'freePar',freePar,...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF',...
    'pullFlag',pullFlag,...
    'SysBudget',24);

[mnu4Sq_tmp90, sin2T4_tmp90] = ...
    KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,0.9,'Mode','New');
mnu4Sq_contour_90 = logspace(log10(min(mnu4Sq_tmp90{1})),log10(max(mnu4Sq_tmp90{1})),5e3);
sin2T4_contour_90 = interp1(mnu4Sq_tmp90{1},sin2T4_tmp90{1},mnu4Sq_contour_90,'lin');

[mnu4Sq_tmp95, sin2T4_tmp95] = ...
    KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,0.95,'Mode','New');
mnu4Sq_contour_95 = logspace(log10(min(mnu4Sq_tmp95{1})),log10(max(mnu4Sq_tmp95{1})),5e3);
sin2T4_contour_95 = interp1(mnu4Sq_tmp95{1},sin2T4_tmp95{1},mnu4Sq_contour_95,'lin');



GetFigure
plot(sin2T4_contour_90,mnu4Sq_contour_90,'LineStyle','-');
hold on;
%plot(sin2T4_tmp90,mnu4Sq_tmp90,'.','MarkerSize',20);
plot(sin2T4_contour_95,mnu4Sq_contour_95,'LineStyle','-.');
set(gca,'YScale','log');
set(gca,'XScale','log');


savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/Files4Thierry/'];
savefile = sprintf('%sSamakContour_%s_%.0feV_%s_%s.mat',savedir,DataType,range,chi2Str,strrep(freePar,' ',''));

if pullFlag<=12
savefile = strrep(savefile,'.mat',sprintf('_pull%.0f.mat',pullFlag));
end
MakeDir(savedir)
save(savefile,'sin2T4_contour_90','mnu4Sq_contour_90','sin2T4_contour_95','mnu4Sq_contour_95');
fprintf('save file to %s \n',savefile);


%%

if contains(freePar,'mNu')
    d = importdata(savefileGrid);
    mNuSq = cell2mat(cellfun(@(x) x.par(1),d.FitResults,'UniformOutput',false));
    GetFigure
    
    mNuSq(abs(mNuSq)>5) = NaN;
    surf(sin2T4,mnu4Sq,mNuSq,'EdgeColor','interp','FaceColor','interp');
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    c =colorbar;
    c.Label.String = sprintf('\\Delta\\chi^2');
    PrettyFigureFormat('FontSize',24)
    view([0 0 1])
    
    c.Label.String = sprintf('{\\itm}_\\nu^2 (eV^2)');
    PrettyFigureFormat('FontSize',24)
    c.Label.FontSize = get(gca,'FontSize')+2;
    grid off
    %c.Limits=[0 zlimMax];
    xlabel('|U_{e4}|^2');
    ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
    zlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
    xlim([1e-03 0.5])
    print(gcf,strrep(savefile,'.mat','.png'),'-dpng','-r500')
end