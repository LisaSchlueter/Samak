
% rebin KNM2TT to match Saenz
FSDdir = [getenv('SamakPath'),'inputs/FSD/'];
SaenzTT = importdata([FSDdir,'FSD_Saenz_T2mod.dat']);
KNM2TT  = importdata([FSDdir,'FSD_KNM2_T2.txt']);

% define energy bins: 
% bin center -> Saenz FSD
% bin edges: use always middle between distance between bin centers as bin edges
Energy    = SaenzTT(:,1);
EnergyL    =  zeros(size(Energy)); % lower edge
EnergyU    =  zeros(size(Energy)); % upper edge

BinWidth   = diff(Energy)./2;
EnergyL(1) =  -1;                           % assume lowest bin covers all smallest excitation energies
EnergyL(2:end) = Energy(2:end)-BinWidth;
EnergyU(1:end-1) = Energy(1:end-1)+BinWidth;
EnergyU(end) = Energy(end)+BinWidth(end);   % assume highted bin covers hightest excitations

SanityPlot = 'OFF';
if strcmp(SanityPlot,'ON')
    % sanity plot
    GetFigure
    plot(Energy,ones(size(Energy)),'x')
    hold on;
    set(gca,'XScale','lin');
    plot(EnergyL,ones(size(Energy)),'o')
    hold on;
    plot(EnergyU,ones(size(Energy)),'.')
    ylim([1-1e-05,1+1e-05]);
end
%%

% loop over energy
ProbKnm2 = zeros(size(Energy));
for i=1:numel(Energy)
%     if i==1
%         SumIdx = find(KNM2TT(:,1)<=Energy(1));
%         ProbKnm2(i) = sum(KNM2TT(SumIdx,2));
%         MeanEnergy(i) = mean(KNM2TT(SumIdx,1));
  %  else      
        SumIdx = find(KNM2TT(:,1) >= EnergyL(i) & KNM2TT(:,1) < EnergyU(i));
        ProbKnm2(i) = sum(KNM2TT(SumIdx,2));
        
 %   end
%    Probtmp = 0;   
%     for j=1:numel(KNM2TT(:,1))
%        if KNM2TT(j,1)>=
%            Probtmp = Probtmp;
%        end
%     end
    
    a=1;
end
%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.6]);
s1 = subplot(4,1,[1:2]);
b1 = bar(Energy,SaenzTT(:,2),'FaceColor',rgb('SkyBlue'));
hold on;
b2 = bar(Energy,ProbKnm2);
set(gca,'XScale','log')
set(gca,'YScale','log')
ylabel('Probability');
xticklabels(''); 
PrettyFigureFormat('FontSize',20);
ax1 = gca;
leg = legend([b2,b1],sprintf('T_2 KNM2'),sprintf('T_2 Saenz2000'),'Location','northwest');
leg.EdgeColor = rgb('Silver');

s2 = subplot(4,1,4);
Diff = ProbKnm2-SaenzTT(:,2);
Mean = 0.5.*(ProbKnm2+SaenzTT(:,2));
RelDiff = Diff./Mean;
p2 = plot(Energy,RelDiff,'.','MarkerSize',10,'Color',rgb('DodgerBlue'));
hold on;
plot(linspace(1e-03,1e3,10),zeros(10,1),'-','Color',rgb('DimGray'),'LineWidth',2);
set(gca,'XScale','log');
ylabel('Rel. diff.');
xlabel('Excitation energy (eV)');
PrettyFigureFormat('FontSize',20);
ax2 = gca;

%%
print('FSD.png','-dpng','-r200');
ax2.Position(4) = 0.2;
ax2.YLabel.Position(1) = 0.0085;
ylim([-0.7,2.3]);
leg = legend(p2,sprintf('(KNM2 - Saenz2000) / mean'),'Location','northwest');
leg.EdgeColor = rgb('Silver');

s3 = subplot(4,1,3);
p3 = plot(Energy,Diff.*1e3,'.','MarkerSize',10,'Color',rgb('ForestGreen'));
hold on;
plot(linspace(1e-03,1e3,10),zeros(10,1),'-','Color',rgb('DimGray'),'LineWidth',2);
linkaxes([s1,s2,s3],'x');
set(gca,'XScale','log');
ylabel(sprintf('Abs. diff.\\times10^3'));
xticklabels('');
PrettyFigureFormat('FontSize',20);
ax3 = gca;
ax3.Position(4) = 0.2;
ax3.YLabel.Position(1) = 0.0085;
xlim([2e-02,300])
ylim([-6 5]);
leg = legend(p3,sprintf('KNM2 - Saenz2000'),'Location','northwest');
leg.EdgeColor = rgb('Silver');

