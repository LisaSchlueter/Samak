% read Stereo data (mat file)
% plot contour
SavePlt = 'ON';
Interp = 'ON';
Type = 'Sensitivity';

   savedir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/'];
   
switch Interp
    case 'OFF'
        if strcmp(Type,'Data')
            matfile = [savedir,'/StereoMaps.mat'];
        else
            matfile = [savedir,'StereoMaps_Sensitivity.mat'];
        end
        
        if exist(matfile,'file')
            load(matfile)
            fprintf('load %s \n',matfile);
        else
            fprintf('file does exist %s \n Run ksn2_ReadStereoChi2Map.m',matfile);
            return
        end
    case 'ON'
        Maxm4Sq = 38^2;
        savefile = sprintf('%sksn2_InterpStereoMini_Max%.0feV2_Min%.2gfeV2.mat',...
            savedir,Maxm4Sq,0.1);
        if strcmp(Type,'Sensitivity')
            savefile = strrep(savefile,'.mat','_Sensitivity.mat');
        end
        d = importdata(savefile);
        fprintf('import %s \n',savefile)
        sin2TSq = d.sin2T4_Osci_cut;
        mNu41Sq = d.mNu4Sq_cut;
        chi2crit = d.chi2Stereocrit_cut;
        chi2  = d.chi2Stereo_cut;
        
end
%% critical delta chi2
GetFigure;
s1 = surf(sin2TSq,mNu41Sq,chi2crit,'EdgeColor','none');%,'EdgeAlpha','interp','EdgeColor','interp');
hold on;
 [C,h]= contour3(sin2TSq,mNu41Sq,chi2crit,[7 8 9],...
        'Color',rgb('White'),'ShowText','on','LineWidth',2.5);
    clabel(C,h,'Color',rgb('White'));
view(2);
set(gca,'Xscale','log');
set(gca,'Yscale','log');
grid off
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('sin(2\\theta)^2'));
ylabel(sprintf('{\\itm}_{41}^2 (eV^2)'));
c = colorbar;
c.Label.String = sprintf('\\Delta\\chi^2_{crit.}');
c.Label.FontSize = get(gca,'FontSize');
t = title('STEREO:  Phys.Rev.D 102 (2020) 052002','FontSize',15,'FontWeight','normal');
pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/plots/'];
if strcmp(SavePlt,'ON')
    MakeDir(pltdir);
    pltname = sprintf('%sStereoChi2crit_Interp%s_%s.png',pltdir,Interp,Type);
    print(gcf,pltname,'-dpng','-r300');
    fprintf('save plot to %s \n',pltname);
end
%%  delta chi2 + contour
chi2Plot = chi2;%NaN.*ones(size(chi2));
chi2Plot(chi2>=chi2crit) = NaN;
GetFigure;
s1 = surf(sin2TSq,mNu41Sq,chi2Plot,'EdgeColor','interp');%,'EdgeAlpha','interp','EdgeColor','interp');
hold on;
%c = contour3(sin2TSq,mNu41Sq,chi2,5.99);
view(2);
grid off
set(gca,'Xscale','log');
set(gca,'Yscale','log');
PrettyFigureFormat('FontSize',22);
xlabel(sprintf('sin(2\\theta)^2'));
ylabel(sprintf('{\\itm}_{41}^2 (eV^2)'));
c = colorbar;
c.Label.String = sprintf('\\Delta\\chi^2');
c.Label.FontSize = get(gca,'FontSize');
if strcmp(Type,'Sensitivity')
    t = title('STEREO Sensitivity:  Phys.Rev.D 102 (2020) 052002','FontSize',15,'FontWeight','normal');
else
    t = title('STEREO:  Phys.Rev.D 102 (2020) 052002','FontSize',15,'FontWeight','normal');
end
if strcmp(SavePlt,'ON')
    MakeDir(pltdir);
    pltname = sprintf('%sStereoChi2_Interp%s_%s.png',pltdir,Interp,Type);
    print(gcf,pltname,'-dpng','-r300');
      fprintf('save plot to %s \n',pltname);
end

