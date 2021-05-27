% read Stereo data (mat file)
% plot contour
SavePlt = 'ON';
Interp = 'ON';

switch Interp
    case 'OFF'
        
        matfile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/StereoMaps.mat'];
        if exist(matfile,'file')
            load(matfile)
            fprintf('load %s \n',matfile);
        else
            fprintf('file does exist %s \n Run ksn2_ReadStereoChi2Map.m',matfile);
            return
        end
    case 'ON'
        
        Maxm4Sq = 40^2;
        DataType = 'Real';
        
        
        savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
        savefile = sprintf('%sksn2_%s_InterpStereo_%s_Max%.0feV2.mat',...
            savedir,DataType,'spline',Maxm4Sq);
        d = importdata(savefile);

        sin2TSq = d.sin2T4_Katrin_Osc;
        mNu41Sq = d.mNu4Sq_Katrin;
        chi2crit = d.chi2critStereo;
        chi2  = d.chi2Stereo;
        
        sin2TSq(~d.InterIdx) = NaN;
        mNu41Sq(~d.InterIdx) = NaN;
        chi2crit(~d.InterIdx) = NaN;
        chi2(~d.InterIdx) = NaN;
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
    print(gcf,[pltdir,sprintf('StereoChi2Crit_Interp%s.png',Interp)],'-dpng','-r300');
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
t = title('STEREO:  Phys.Rev.D 102 (2020) 052002','FontSize',15,'FontWeight','normal');
if strcmp(SavePlt,'ON')
    MakeDir(pltdir);
    print(gcf,[pltdir,sprintf('StereoChi2_Interp%s.png',Interp)],'-dpng','-r300');
end

