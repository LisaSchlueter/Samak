maxE = 500;
ELossBinStep = 0.04;
minE=-maxE; NbinE = (maxE-minE)/ELossBinStep;
E = minE:ELossBinStep:maxE;
Estep = E(2) - E(1);
%E = linspace(minE,maxE,NbinE); Estep = E(2) - E(1);
savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
savename = sprintf('%sELoss_KatrinT2_%.2feVbinning.mat',savedir,ELossBinStep);
if exist(savename,'file') 
    load(savename);
else
    A = ref_FakeRun_KNM2_RFcomparison;
    [~, ElossFunctions] = A.ComputeELossFunction('E',E);
    
    Esave = E(E>=0 & E<90);
    ElossFunctionssave       = ElossFunctions(:,E>=0 & E<90);
    ElossFunctionssave(2,:) = ElossFunctionssave(2,:).*Estep;
    ElossFunctionssave(3,:) = ElossFunctionssave(3,:).*Estep.^2;
    ElossFunctionssave(4,:) = ElossFunctionssave(4,:).*Estep.^3;
    ElossFunctionssave(5,:) = ElossFunctionssave(5,:).*Estep.^4;
    ElossFunctionssave(6,:) = ElossFunctionssave(6,:).*Estep.^5;
    ElossFunctionssave(7,:) = ElossFunctionssave(7,:).*Estep.^5;
    save(savename,'ElossFunctions','E','Estep');
    Write2Txt('filename',strrep(savename,'.mat',''),... % txt file
        'variable',[Esave',ElossFunctionssave']',...
        'nCol',8);
end

%% display
f1 = ElossFunctions(1,:);            % one scattering
f2 = ElossFunctions(2,:).*Estep;     % two scatterings
f3 = ElossFunctions(3,:).*Estep.^2;  % three scatterings
f4 = ElossFunctions(4,:).*Estep.^3;  % fours scatterings
f5 = ElossFunctions(5,:).*Estep.^4;  % five scatterings


fig5 = figure('Renderer','opengl');
set(fig5, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
plot(E,f1,'Linewidth',3,'Color',rgb('FireBrick'));
hold on;
plot(E,f2,'Linewidth',3,'Color',rgb('Orange'),'LineStyle',':');
plot(E,f3,'Linewidth',3,'Color',rgb('RoyalBlue'),'LineStyle','--');
plot(E,f4,'Linewidth',3,'Color',rgb('CadetBlue'),'LineStyle','-.');%plot(E,f5(E).*EBinSize,'Linewidth',3,'Color',rgb('Black'),'LineStyle','-.');
plot(E,f5,'Linewidth',3,'Color',rgb('SlateGray'),'LineStyle',':');

PrettyFigureFormat;
set(gca,'FontSize',24);
xlim([0 90]);
ylim([0 0.27])
xlabel(['energy loss ',char(949),' (eV)']);
ylabel(['probability f(',char(949),')']);
leg = legend(string(1:4)); legend boxoff;
leg.FontSize = 24;
leg.Title.String = 'number of scatterings';
% save_name = sprintf('../../studies/ft_CovarianceMatrices/plots/EnergyLossFunctions/ELossFunctions_%s',ELossFlag);
% print(fig5,[save_name,'.png'],'-dpng','-r300');
% publish_figurePDF(fig5,[save_name,'.pdf']);

