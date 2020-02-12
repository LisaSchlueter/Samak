% T-minus FSD
% script for display 
%%  preliminary T-minos ion final state distribution
Tm = importdata('FSD_Tminus_unshifted.mat');
E    = Tm(:,1); % excitation energy (unshifted) from https://www.researchgate.net/publication/252490208_Final-State_Spectrum_of_3He_after_beta-_Decay_of_Tritium_Anions_T-
Prob = Tm(:,2); % probability
BSIndex = 15;  % bound state threshold (both electrons bound)

%% shift FSD to correct energy reference
EShift = -12.17-min(E); %shift ground state by difference of endpoint (T^2 vs T_2)
FSD = [E+EShift,Prob];
save([getenv('SamakPath'),'inputs/FSD/FSD_Tminus.mat'],'FSD');
%% plot
PlotMode = 'Shift'; % show shited or unshifted fsd

figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
switch PlotMode
    case 'NoShift'
        pES = plot(E(BSIndex+1:end),Prob(BSIndex+1:end)*100,'-','LineWidth',3,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
        hold on;
        pBS = plot(E(1:BSIndex),Prob(1:BSIndex)*100,'o','LineWidth',3,'MarkerFaceColor',rgb('DodgerBlue'),'Color',rgb('DodgerBlue'));
    case 'Shift'
        pES = plot(FSD(BSIndex+1:end,1),FSD(BSIndex+1:end,2)*100,'-','LineWidth',3,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
        hold on
        pBS = plot(FSD(1:BSIndex,1),FSD(1:BSIndex,2)*100,'o','LineWidth',3,'MarkerFaceColor',rgb('DodgerBlue'),'Color',rgb('DodgerBlue'));
end
hold off
PrettyFigureFormat('FontSize',24);
xlabel('Excitation energy (eV)');
ylabel('Probability (%)');
xlim([-20,150])
set(gca,'YScale','log');
leg = legend([pBS,pES],'Bound states','Continous states');
leg.Title.String = sprintf('FSD (T^-\\rightarrow^3He)'); 
leg.FontSize = get(gca,'FontSize');
legend boxoff

%% save plot
savefile = [getenv('SamakPath'),sprintf('knm2ana/knm2_Plasma/plots/TminusFSD_%s.pdf',PlotMode)];
%export_fig(gcf,savefile,'-painters');
%pfprintf('save plot to %s \n',savefile)
%print(gcf,savefile,'-r400','-dpng');