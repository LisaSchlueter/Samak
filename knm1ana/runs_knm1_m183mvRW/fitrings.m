% Fit all 175mV Runs (up to 9/04/19) Ringwise
%% Init
MyRings=1:12;
RunList = 'KNM1_m149mvRW';
RWlabel = strrep(strrep(RunList,'KNM1_',''),'vRW','V');
recomputeFlag = 'ON';
range = 90;
filepath = sprintf('./results/');
if ~exist(filepath,'dir')
    system(['mkdir ',filepath]);
end
filename = sprintf('FitRings%s_%.0frange.mat',RWlabel,range);
if exist([filepath,filename],'file') && strcmp(recomputeFlag,'OFF')
    load([filepath,filename]);
else
%%
if range==90
    exclDataStart = 2;
elseif range==40
    exclDataStart=12;
end
Rings = MyRings;
nRings =numel(Rings);
E0 = zeros(nRings,1);
E0Err = zeros(nRings,1);
chi2min = zeros(nRings,1);
dof   = zeros(nRings,1);
M = arrayfun(@(x) MultiRunAnalysis('RunList',RunList,'DataType','Real',...
'AnaFlag','Ring','RingList',x,'exclDataStart',exclDataStart),Rings);
%% spectra
TBDISrings  = cell2mat(arrayfun(@(x) x.RunData.TBDIS,M,'UniformOutput',0));
qUrings     = cell2mat(arrayfun(@(x) x.RunData.qU,M,'UniformOutput',0));
qUfracrings = cell2mat(arrayfun(@(x) x.RunData.qUfrac,M,'UniformOutput',0));
Timerings   = cell2mat(arrayfun(@(x) x.RunData.TimeSec,M,'UniformOutput',0));
StackedRuns = cell2mat(arrayfun(@(x) numel(x.StackedRuns),M,'UniformOutput',0));
RateRings = TBDISrings./(qUfracrings.*Timerings);
nPix  = cell2mat(arrayfun(@(x) numel(x.PixList),M,'UniformOutput',0));

%% Fit
parfor i=1:nRings
    M(i).Fit;
    E0(i) = M(i).FitResult.par(2);
    E0Err(i) = M(i).FitResult.err(2);
    chi2min(i) = M(i).FitResult.chi2min;
    dof(i) = M(i).FitResult.dof;
end
%% save
save([filepath,filename],...
    'Rings','nRings','E0','E0Err','chi2min','M',...
    'TBDISrings','qUrings','qUfracrings','range','RateRings','Timerings','exclDataStart','StackedRuns','nPix');
end
%% linear fit
meanE0 =wmean(E0,1./E0Err.^2);
[par, err, chi2min_linfit,dof_linFit] =linFit(Rings',E0-meanE0,E0Err);

%% plot fit result
fig2 = figure('Renderer','opengl');
set(fig2,'units','normalized','pos',[0.1, 0.1,0.7,0.7]);
%b = bar(Rings,E0-meanE0,'FaceColor',rgb('SkyBlue'),'LineStyle','none');
e = errorbar(Rings,E0-meanE0,E0Err,'o',...
    'LineWidth',2,'LineStyle','none','MarkerSize',8,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('SlateGray'));
hold on;
plot(linspace(0.5,12.5,nRings),zeros(nRings,1),'--','Color',rgb('SlateGray'),'LineWidth',2);
l = plot(Rings,par(1).*Rings+par(2),'-','Color',rgb('SkyBlue'),'LineWidth',3);
hold off;
PrettyFigureFormat;
ylabel('E_0 - <E_0> (eV)')
xlabel('ring')
set(gca,'FontSize',20);
title(sprintf('%.0f stacked runs (%.0f - %.0f) - %.0f pixels total',numel(M(1).StackedRuns),M(1).RunList(1),M(1).RunList(end),118));
leg = legend([e,l],...
    sprintf('%.0f eV range',range),...
    sprintf('linear fit slope: (%.1f \\pm %.1f) meV @ \\chi2 = %.1f / %.0f dof',par(1)*1e3,err(1)*1e3,chi2min_linfit,dof_linFit));
legend boxoff
xlim([0.5,12.5])
leg.Location = 'northwest';
if ~exist('./plots/','dir')
    system('mkdir ./plots/');
end

savename = sprintf('./plots/E0ringwise_%sRuns_%.0frange.png',RWlabel,range);
print(savename,'-dpng','-r500');
%%
%Write2Txt('filename',[filepath,strrep(filename,'.mat','')],'variable',[E0-meanE0,E0Err]','variableName','E0-meanE0      E0Err','nCol',2);
%% plot background
fig3 = figure('Renderer','opengl');
set(fig3,'units','normalized','pos',[0.1, 0.1,0.7,0.7]);
BkgqU = qUrings(:,1)>18575;
BKG   = RateRings(BkgqU,:);
BKGRate = BKG./(qUfracrings(BkgqU,:).*Timerings);
%plot(Rings,BKG./nPix*1e3,'o-');
%hold on;
plot(Rings,mean(BKGRate./nPix)*1e3,'o-',...
    'LineWidth',2,'Color',rgb('SlateGray'),'MarkerFaceColor',rgb('SlateGray'),'MarkerSize',8);
ylabel('rate per npixels (cps)');
hold on;
scaling = 0.0004;
plot(Rings,scaling.*mean(RateRings)./nPix,'o-',...
    'LineWidth',2,'Color',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'MarkerSize',8);
leg = legend('mean background (5 subruns)',sprintf('mean rate \\cdot %g (all subruns)',scaling)); legend boxoff
xlabel('ring');
PrettyFigureFormat;
set(gca,'FontSize',20);
xlim([0.5,12.5])

print(sprintf('./plots/BackgroundRateRings_%sRuns_%.0feVrange.png',RWlabel,range),'-dpng','-r500');
%% plot spectrum
% Mode = 'Counts';
% fig1 = figure('Renderer','opengl');
% set(fig1,'units','normalized','pos',[0.1, 0.1,0.7,0.7]);
% 
% switch Mode
%     case 'Counts'
%         p = plot(qUrings-18575,TBDISrings,'o-',...
%         'LineWidth',2,'MarkerSize',5);
%         ylabel('counts');
%     case 'Rate'
%         p = plot(qUrings-18575,RateRings,'o-','LineWidth',2,'MarkerSize',5);
%         ylabel('count rate (cps)');
% end
% 
% set(p,{'color'},num2cell(jet(12),2));
% PrettyFigureFormat;
% leg = legend(string(Rings));
% legend boxoff
% leg.Title.String = 'ring';
% xlabel('retaring potential - 18575(eV)');
% title(sprintf('%.0f stacked runs (%.0f - %.0f)',numel(M(1).StackedRuns),M(1).RunList(1),M(1).RunList(end)));
% xlim([-94 47]);
% set(gca,'FontSize',20);
%savename = sprintf('../KNM1_RingAnalysis/plots/TBDIS_%sRuns.png',RWlabel)
% print(savename,'-dpng','-r500');

%% FPD Viewer
% Fill Results in Pixel Array
% E0pixel=ones(1,148)*NaN;
% for r=2:11
% E0pixel([M(r).ModelObj.ring{r}]) = [E0(r)-meanE0];
% end
% FPDViewer(E0pixel,'ReDrawSkeleton','OFF')
% set(gca,'visible','off')
% colormap(cool)