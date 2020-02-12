% script in plot neutrino mass shift and endpoint shift for different twins

RunList = 'KNM1';
exclDataStart_all = [17,14];
ReFit = 'ON';
nTwins = 7;
mNuSq    = zeros(numel(exclDataStart_all),nTwins); mNuSqErr = zeros(numel(exclDataStart_all),nTwins);
E0       = zeros(numel(exclDataStart_all),nTwins); E0Err    = zeros(numel(exclDataStart_all),nTwins);
chi2min  = zeros(numel(exclDataStart_all),nTwins);
dof = zeros(numel(exclDataStart_all),1);

for i=1:numel(exclDataStart_all)
[mNuSq(i,:),mNuSqErr(i,:),E0(i,:),E0Err(i,:),~,~,~,~,dof(i),chi2min(i,:)] = FitTwins('RunList',RunList,...
    'ReFit',ReFit,'exclDataStart',exclDataStart_all(i),'SavePlot','OFF');
end

%%
Parameter = 1;
switch Parameter
    case 1 % nu-mass
        y = mNuSq;
        yErr = mNuSqErr;
        ystr = sprintf('m^2_\\nu shift (eV^2)');
    case 2 % endpoint
        y    = E0;
        yErr = E0Err;
        ystr = sprintf('E_0 shift (eV)');
    case 3 % chi2
        y = chi2min;%1-chi2cdf(chi2min,repmat(dof,1,nTwins));
        yErr = zeros(numel(exclDataStart_all),nTwins);
        ystr = sprintf('\\chi2');
end

twinLabel = {'twin',sprintf('same \\rhod'),'same T_2','same qUf','same qU','same qU qUf','all same'};%,'all same -qUf'
PlotColor   = {rgb('IndianRed');rgb('CadetBlue');rgb('DarkGoldenRod')};
MarkerColor = {rgb('FireBrick');rgb('RoyalBlue');rgb('Sienna')};
lineArg = {'o-','LineWidth',4,'MarkerSize',10};
x = 1:numel(twinLabel);

% plot
fig1 = figure(1);
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
legstr = cell(numel(exclDataStart_all),1);

for i=1:numel(exclDataStart_all)
plot(x,y(i,:),lineArg{:},'MarkerFaceColor',MarkerColor{i},'Color',PlotColor{i});
hold on;
legstr{i} = sprintf('%.0f eV',Convert2Range(exclDataStart_all(i)));
end
hold off;

leg = legend(legstr{:}); legend boxoff
leg.Title.String = 'range'; leg.FontSize = 22;
leg.Location = 'northwest';
xticks(1:numel(x));
xticklabels(twinLabel);
xtickangle(35)
yl = ylabel(ystr); yl.FontSize = 22;
PrettyFigureFormat;
grid on;
% save
savedir = [getenv('SamakPath'),'knm1ana/knm1_Stacking/plots/'];
if ~exist(savedir,'dir')
    system(['mkdir -p ',savedir]);
end

print(gcf,[savedir,sprintf('Fit2Twins_%.0f.png',Parameter)],'-dpng','-r450');
fprintf('plot saved to %s \n',[savedir,sprintf('Fit2Twins_%.0f.png',Parameter)]);

function range = Convert2Range(excldatastart)
switch excldatastart
    case 2
        range = 90;
    case 14
        range = 40;
    case 17
        range = 30;
end
end