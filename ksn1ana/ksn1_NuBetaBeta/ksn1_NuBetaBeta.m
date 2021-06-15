SanityPlot = 'ON';

%% load files
range = 40;
DataType = 'Real';
freePar = 'mNu E0 Norm Bkg';
pullFlag = 15;
switch pullFlag
    case 99
        pullStr = '';
    case 15
        pullStr = sprintf('_mNuSqPull1eV2');
end
savedir = [getenv('Samak3.0'),'ksn1ana/ksn1_NuBetaBeta/results/'];

file1 = sprintf('%sContour_Real_%s_40eV.mat',savedir,'E0NormBkg');
file2 = sprintf('%sContour_Real_%s_40eV.mat',savedir,'mNuE0NormBkg');
file3 = sprintf('%sContour_Real_%s_40eV_mNuSqPull1eV2.mat',savedir,'mNuE0NormBkg');

d1 = importdata(file1);
d2 = importdata(file2);
d3 = importdata(file3);
%% sanity plot (contours)
if strcmp(SanityPlot,'ON')
GetFigure;
plot(d1.sinT4Sq,d1.mNu4Sq,'LineWidth',1.5);
hold on;
plot(d2.sinT4Sq,d2.mNu4Sq,'LineWidth',1.5);
plot(d3.sinT4Sq,d3.mNu4Sq,'LineWidth',1.5);
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',20);
leg = legend(sprintf('i) {\\itm}_\\nu^2 fix'),sprintf('ii) {\\itm}_\\nu^2 free'),sprintf('iii) {\\itm}_\\nu^2 free, \\sigma({\\itm}_\\nu^2) = 1 eV^2'));
legend boxoff
leg.Location = 'southwest';
ylim([1 2e3]); 
xlim([1e-02 0.5]);
hold off;
end

%%
GetFigure;
h1 = histogram(d1.sinT4Sq.*sqrt(d1.mNu4Sq),'Normalization','probability',...
    'EdgeColor',rgb('Orange'),'FaceColor',rgb('Orange'),'FaceAlpha',0.2,'LineWidth',1.5);
hold on;
h2 = histogram(d2.sinT4Sq.*sqrt(d2.mNu4Sq),'Normalization','probability',...
    'EdgeColor',rgb('SteelBlue'),'FaceColor',rgb('SkyBlue'),'FaceAlpha',0.2,'LineWidth',1.5,'LineStyle',':');
ylimits = ylim;
%p1 = plot(d1.sinT4Sq_bf.*d1.mNu4Sq_bf.*ones(5,1),linspace(1e-03,1,5),'-','LineWidth',2,'Color',h1.EdgeColor);
%p2 = plot(d2.sinT4Sq_bf.*d2.mNu4Sq_bf.*ones(5,1),linspace(1e-03,1,5),'-','LineWidth',2,'Color',h2.EdgeColor);

set(gca,'YScale','log');
ylabel('Frequency');
ylim([1e-03,1]);
xlabel(sprintf('{\\itm}_4 \\times |{\\itU}_{e4}|^2 (eV)'))
PrettyFigureFormat('FontSize',20);
leg = legend(sprintf('i) {\\itm}_\\nu^2 fix    (best-fit: {\\itm}_4 \\times |{\\itU}_{e4}|^2 = %.1f eV)',d1.sinT4Sq_bf.*sqrt(d1.mNu4Sq_bf)),...
    sprintf('ii) {\\itm}_\\nu^2 free (best-fit: {\\itm}_4 \\times |{\\itU}_{e4}|^2 = %.1f eV)',d2.sinT4Sq_bf.*sqrt(d2.mNu4Sq_bf)));
leg.Title.String = sprintf('Exclusion at 95%% C.L.');
leg.Title.FontWeight = 'normal';
legend boxoff
leg.Location = 'northeast';

plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir)
plotname = sprintf('%sNuBetaBeta_Product_m4sinT.png',plotdir);
print(plotname,'-dpng','-r350');

%%
 close all
% 
% surf(d1.sinT4Sq,d1.mNu4Sq,d1.sinT4Sq.*d1.mNu4Sq','EdgeColor','none')
% hold on;
% plot3(d1.sinT4Sq,d1.mNu4Sq,1e2*ones(numel(d1.mNu4Sq),1),...
%                     '-','MarkerSize',9,'Color',rgb('White'),'LineWidth',3);

                
nGrid = numel(d1.mNu4Sq);
m4  = d1.mNu4Sq;%logspace(log(min(d1.mNu4Sq)),log(max(d1.mNu4Sq)),nGrid);
sin = d1.sinT4Sq;%logspace(log(min(d1.sinT4Sq)),log(max(d1.sinT4Sq)),nGrid);

%logicsin = logical(sin<d1.sinT4Sq');%sin<d1.sinT4Sq;
%logicTot = ~logical(logicsin);%logical(sin<d1.sinT4Sq' & m4'<d1.mNu4Sq) ;%~logical(logicm4.*logicsin);

z = sqrt(m4)'.*sin;
z(sin>d1.sinT4Sq') = NaN;
zCutOff = 1.2;
z(z>zCutOff) = zCutOff;

%a = z;
%a(a>=0.199 & a<=0.201) = 99;

surf(sin',m4,z,'EdgeColor','none');
hold on;
%surf(sin',m4,a,'EdgeColor','White');
surf(sin',m4,-1.*ones(nGrid),'EdgeColor','none','FaceColor',rgb('LightGray'));
plot3(d1.sinT4Sq,d1.mNu4Sq,1e3*ones(numel(d1.mNu4Sq),1),...
                    '-','MarkerSize',9,'Color',rgb('White'),'LineWidth',2);
                x = logspace(min(log(d1.sinT4Sq)),log(0.5),10);
pt = plot3(x,(0.2./x).^2,1e3*ones(10,1),...
'--','Color',rgb('Black'),'LineWidth',1.5);
t = text(0.028,35,sprintf('{\\itm}_4 \\times |{\\itU}_{e4}|^2 = 0.2 eV'),'Rotation',-42.5,'FontSize',15);
view(2);
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',20);
grid off;
xlim([min(sin) max(sin)]); 
ylim([min(m4) 1e3]);%max(m4)]);
c = colorbar;
c.Limits = [0 zCutOff];
c.Label.String = sprintf('{\\itm}_4 \\times |{\\itU}_{e4}|^2 (eV)');
c.Label.FontSize = 19;
plotname = sprintf('%sNuBetaBeta_ProductMap_m4sinT.png',plotdir);
print(plotname,'-dpng','-r350');