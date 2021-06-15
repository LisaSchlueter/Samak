InterpMode = 'spline';
Maxm4Sq = 40^2;%40^2;
DataType = 'Real';
 SanityPlt = 'ON';

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
savefile = sprintf('%sksn2_%s_InterpStereo_%s_Max%.0feV2.mat',...
    savedir,DataType,InterpMode,Maxm4Sq);

if exist(savefile,'file')
    fprintf('load file %s \n',savefile);
    load(savefile)
else
    fprintf('file not found. run: ksn2_InterpStereo.m \n')
    return
end

if strcmp(SanityPlt,'OFF')
    figure('Units','normalized','Position',[0.1,0.1,0.8,0.8]);
    % katrin + stereo
    s1=  subplot(2,2,1);
    [pcombi] = surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2Combi,'EdgeColor','none');
    hold on;
    pcontour =plot3(sin2T4_Katrin_contour,mNu4Sq_Katrin_contour,ones(numel(mNu4Sq_Katrin_contour),1).*1e4,...
        'LineWidth',2,'Color',rgb('White'));
    pbf = plot3(sin2T4_Katrin_bf,mNu4Sq_Katrin_bf,1e4,'wx','LineWidth',2);
    grid off
    view(2)
    set(gca,'Xscale','log')
    set(gca,'yscale','log')
    xlabel(sprintf('|{\\itU}_{e4}|^2'));
    ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
    PrettyFigureFormat;
    c1 = colorbar;
    c1.Label.String = sprintf('\\Delta\\chi^2');
    c1.Label.FontSize = get(gca,'FontSize');
    %title('KATRIN','FontSize',15,'FontWeight','normal');
    leg = legend([pcombi,pcontour],'KATRIN + STEREO','KATRIN only exclusion at 95% C.L.','Location','southwest');
    PrettyLegendFormat(leg);
    %% stereo only
    s2 = subplot(2,2,2);
   chi2StereoPlt = NaN.*zeros(size(chi2_Katrin));
   chi2StereoPlt(logical(InterIdx)) = chi2Stereo(logical(InterIdx));  
    surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2StereoPlt,'EdgeColor','none') %
    grid off
    view(2)
    set(gca,'Xscale','log')
    set(gca,'yscale','log')
    xlabel(sprintf('|{\\itU}_{e4}|^2'));
    ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
    PrettyFigureFormat;
    c2 = colorbar;
    c2.Label.String = sprintf('\\Delta\\chi^2');
    c2.Label.FontSize = get(gca,'FontSize');
    leg = legend(sprintf('STEREO'),'Location','southwest'); %(\\chi^2_{crit.} map 95%% C.L.)
    PrettyLegendFormat(leg);
    
     
    linkaxes([s1,s2],'xy');
    ylim([0.1 40^2]);
    xlim([1e-03 0.5])
    %% katrin + stereo: zoom at stereo region
    s3 = subplot(2,2,3);
    
    tmp = chi2Combi;
    chi2combiPlt  = NaN.*zeros(size(chi2_Katrin));
    chi2combiPlt(logical(InterIdx)) = tmp(logical(InterIdx));
    
    surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2combiPlt,'EdgeColor','none') %
    hold on;
    [C,h]= contour3(sin2T4_Katrin,mNu4Sq_Katrin,chi2combiPlt,[5, 7 10 ],...
        'Color',rgb('White'),'ShowText','on');
    clabel(C,h,'Color',rgb('White'));
      grid off
    view(2)
    set(gca,'Xscale','log')
    set(gca,'yscale','log')
    xlabel(sprintf('|{\\itU}_{e4}|^2'));
    ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
    PrettyFigureFormat;
    leg = legend('KATRIN + STEREO (Zoom)','Location','southwest');
    PrettyLegendFormat(leg);
    c2 = colorbar;
    c2.Label.String = sprintf('\\Delta\\chi^2');
    c2.Label.FontSize = get(gca,'FontSize');
    
    %% stereo chi2 crit map
    s4 = subplot(2,2,4);
    surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2critStereo,'EdgeColor','none') %
    hold on;
    [C,h]= contour3(sin2T4_Katrin,mNu4Sq_Katrin,chi2critStereo,[5 7 8 ],...
        'Color',rgb('White'),'ShowText','on');
    clabel(C,h,'Color',rgb('White'));
    grid off
    view(2)
    set(gca,'Xscale','log')
    set(gca,'yscale','log')
    xlabel(sprintf('|{\\itU}_{e4}|^2'));
    ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
    PrettyFigureFormat;
    leg = legend(sprintf('KATRIN (WILKS) + STEREO (\\chi^2_{crit} map)'),'Location','southwest');
    PrettyLegendFormat(leg);
    c2 = colorbar;
    c2.Label.String = sprintf('\\Delta\\chi^2_{crit.}');
    c2.Label.FontSize = get(gca,'FontSize');
    
    linkaxes([s3,s4],'xy');
    ylim([0.1 11]);
    xlim([0.002 0.5])
    
    pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/plots/'];
    MakeDir(pltdir);
    print(gcf,[pltdir,'KATRINandStereo.png'],'-dpng','-r300');
end

%% KATRIN + STEREO: critical chi^2
InterIdx = logical(InterIdx);
GetFigure%chi2_Katrin+chi2Stereo
surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2critStereo,'EdgeColor','none');
% hold on;
% pcontour =plot3(sin2T4_Katrin_contour,mNu4Sq_Katrin_contour,ones(numel(mNu4Sq_Katrin_contour),1).*1e4,...
%     'LineWidth',2,'Color',rgb('White'));
% pbf = plot3(sin2T4_Katrin_bf,mNu4Sq_Katrin_bf,1e4,'wx','LineWidth',2);
view(2)
grid off;
set(gca,'Xscale','log')
set(gca,'yscale','log')
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
PrettyFigureFormat;
c2 = colorbar;
c2.Label.String = sprintf('\\Delta\\chi_{crit.}^2');
c2.Label.FontSize = get(gca,'FontSize');
ylim([0 40^2]);
xlim([1e-03 0.5])

 pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/plots/'];
    MakeDir(pltdir);
    print(gcf,[pltdir,'KATRINandStereoCrit_Combi.png'],'-dpng','-r300');
%% KATRIN + STEREO: chi^2
    InterIdx = logical(InterIdx);
    GetFigure%chi2_Katrin+chi2Stereo
    surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2Combi,'EdgeColor','none');
    hold on;
    pcontour =plot3(sin2T4_Katrin_contour,mNu4Sq_Katrin_contour,ones(numel(mNu4Sq_Katrin_contour),1).*1e4,...
        'LineWidth',2,'Color',rgb('White'));
    pbf = plot3(sin2T4_Katrin_bf,mNu4Sq_Katrin_bf,1e4,'wx','LineWidth',2);
    view(2)
    grid off;
    set(gca,'Xscale','log')
    set(gca,'yscale','log')
    xlabel(sprintf('|{\\itU}_{e4}|^2'));
    ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
    PrettyFigureFormat;
    c2 = colorbar;
    c2.Label.String = sprintf('\\Delta\\chi^2');
    c2.Label.FontSize = get(gca,'FontSize');
    ylim([0 40^2]);
    xlim([1e-03 0.5])
    
    pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/plots/'];
    MakeDir(pltdir);
    print(gcf,[pltdir,'KATRINandStereo_Combi.png'],'-dpng','-r300');
%% combi exclusion
% best fit (in shared parameter space)
[row,col]= find(chi2Combi(InterIdx)==min(min(chi2Combi(InterIdx))));
sin2T4_bf_tmp = sin2T4_Katrin(InterIdx);
sin2T4_combi_bf =sin2T4_bf_tmp(row,col);
mNu4Sq_tmp = mNu4Sq_Katrin(InterIdx);
mNu4Sq_combi_bf = mNu4Sq_tmp(row,col);

GetFigure;

chi2Combi(chi2Combi>chi2critStereo) = NaN;
[pcombi] = surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2Combi,'EdgeColor','none');
hold on;
pbf = plot3(sin2T4_Katrin_bf,mNu4Sq_Katrin_bf,1e4,'wx','LineWidth',2);
pbf_combi = plot3(sin2T4_combi_bf,mNu4Sq_combi_bf,1e4,'ko','LineWidth',2);
hold on;
grid off
view(2)
set(gca,'Xscale','log')
set(gca,'yscale','log')
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
PrettyFigureFormat;
c1 = colorbar;
c1.Label.String = sprintf('\\Delta\\chi^2');
c1.Label.FontSize = get(gca,'FontSize');
%title('KATRIN','FontSize',15,'FontWeight','normal');
leg = legend([pcombi,pbf,pbf_combi],...
    'KATRIN + STEREO at 95% C.L.','Best fit KATRIN','Best fit KATRIN + STEREO','Location','southwest');
PrettyLegendFormat(leg);
ylim([0.1 40^2])
xlim([1e-03 0.5])
pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/plots/'];
MakeDir(pltdir);
print(gcf,[pltdir,'KATRINandStereo_Combi.png'],'-dpng','-r300');

%% combi exclusion zoom
% best fit
[row,col]= find(chi2Combi==min(min(chi2Combi)));
sin2T4_bf = sin2T4_Katrin(row,col);
mNu4Sq_bf = mNu4Sq_Katrin(row,col);

GetFigure;

%Combi  =chi2_Katrin+chi2Stereo;
chi2Combi(~InterIdx)=NaN;
%Combi(Combi>100) = NaN;
[pcombi] = surf(sin2T4_Katrin,mNu4Sq_Katrin,chi2Combi,'EdgeColor','none');
 hold on;
    [C,h]= contour3(sin2T4_Katrin,mNu4Sq_Katrin,chi2Combi,[12  100],...
        'Color',rgb('Silver'),'ShowText','on','LineWidth',1.5);
    clabel(C,h,'Color',rgb('Silver'),'FontSize',12);
        [C,h]= contour3(sin2T4_Katrin,mNu4Sq_Katrin,chi2Combi,[10  10],...
        'Color',rgb('HotPink'),'ShowText','on','LineWidth',2);
    clabel(C,h,'Color',rgb('HotPink'),'FontWeight','bold','FontSize',12);
    
hold on;
grid off
view(2)
set(gca,'Xscale','log')
set(gca,'yscale','log')
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('\\Deltam_{4}^2 (eV^2)'));
PrettyFigureFormat;
c1 = colorbar;
c1.Label.String = sprintf('\\Delta\\chi^2');
c1.Label.FontSize = get(gca,'FontSize');
leg = legend(pcombi,'KATRIN + STEREO','Location','southwest');
PrettyLegendFormat(leg);
ylim([0.1 10])
xlim([0.0025 0.5])
 pltdir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/plots/'];
 MakeDir(pltdir);
print(gcf,[pltdir,'KATRINandStereo_CombiZoom.png'],'-dpng','-r300');