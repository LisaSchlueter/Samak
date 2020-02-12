clear all;
addpath(genpath('../../../Samak2.0'));

A=InitKATRIN('PS_Wein93','OFF','TD','Flat20','BKG_RateAllFPDSec',1e-3);
B1=InitKATRIN('PS_Wein93','ON','TD','Flat20','BKG_RateAllFPDSec',1e-3,'PS_Wein93Coeff',0.716);
B2=InitKATRIN('PS_Wein93','ON','TD','Flat20','BKG_RateAllFPDSec',1e-3,'PS_Wein93Coeff',0.74);
B3=InitKATRIN('PS_Wein93','ON','TD','Flat20','BKG_RateAllFPDSec',1e-3,'PS_Wein93Coeff',0.76);

%%
figure(1)

subplot(2,1,1)
m=-1;
A.mnuSq_i=m; A.ComputeTBDDS();  
h1 = plot(-A.Q+A.Te,A.PhaseSpace,...
    'LineWidth',2,'LineStyle','-','Color','Blue');
hold on
B1.mnuSq_i=m; B1.ComputeTBDDS();  
h2 = plot(-A.Q+B1.Te,B1.PhaseSpace,...
    'LineWidth',2,'LineStyle','-','Color','Red');
B2.mnuSq_i=m; B2.ComputeTBDDS();  
h3 = plot(-A.Q+B2.Te,B2.PhaseSpace,...
    'LineWidth',2,'LineStyle','-','Color','Black');
B3.mnuSq_i=m; B3.ComputeTBDDS();  
h4 = plot(-A.Q+B3.Te,B3.PhaseSpace,...
    'LineWidth',2,'LineStyle','-','Color','Magenta');
hold off
grid on
xlabel('Te - E_0 (eV)','FontSize',12);
ylabel('Phase Space (a.u.)','FontSize',12);
set(gca,'FontSize',12);
set(gca,'yscale','log');
%axis([-2 0.5 2e9 2e11]);
a = legend([h1 h2 h3 h4],...
    'Regular - m^2=-1 eV^2',...
    'Modified (0.716) - m^2=-1 eV^2','Modified (0.73)- m^2=-1 eV^2','Modified (0.76) - m^2=-1 eV^2',...
    'Location','SouthWest') ; %legend(a,'boxoff','FontSize',10);
%PrettyFigureFormat

subplot(2,1,2)
m=0;
A.mnuSq_i=m; A.ComputeTBDDS();  
h1 = plot(-A.Q+A.Te,A.PhaseSpace,...
    'LineWidth',2,'LineStyle','-','Color','Blue');
hold on
B1.mnuSq_i=m; B1.ComputeTBDDS();  
h2 = plot(-A.Q+B1.Te,B1.PhaseSpace,...
    'LineWidth',2,'LineStyle','-','Color','Red');
hold off
grid on
xlabel('Te - E_0(eV)','FontSize',12);
ylabel('Phase Space (a.u.)','FontSize',12);
set(gca,'FontSize',12);
set(gca,'yscale','log');
%axis([-2 0.5 2e9 2e11]);
    a = legend([h1 h2],'Regular - m^2=0','Modified - m^2=0','Location','NorthEast') ;% legend(a,'boxoff','FontSize',10);
%PrettyFigureFormat

mtit(sprintf('KATRIN - %s - %g y - %g mcps',A.TD,A.TimeYear,A.BKG_RateSec_i));

