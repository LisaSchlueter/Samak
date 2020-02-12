%
% Unveiling the FSD Onset
%
% Thierry Lasserre
% April 2017
%

addpath('./tools','-begin');

TimeYear = 10/365.25;

% Without FSD
clear Foff; global Foff ; 
Foff=InitKATRIN('TimeYear',TimeYear,'TTFSD','OFF','DTFSD','OFF','HTFSD','OFF');
Foff.ComputeTBDDS(); Foff.ComputeTBDIS(); 
%AddStatFluctTBDIS(Foff);
Foff.DisplayTDBInfo();

% With FSD - With Electronics Excited States
clear Fon; global Fon ; 
Fon=InitKATRIN('TimeYear',TimeYear,'TTFSD','SAENZ','DTFSD','OFF','HTFSD','OFF');
Fon.ComputeTBDDS(); Fon.ComputeTBDIS(); 
%AddStatFluctTBDIS(Fon);
Fon.DisplayTDBInfo();

% With FSD - Without Electronics Excited States
clear Fonnoee; global Fonnoee ; 
Fonnoee=InitKATRIN('TimeYear',TimeYear,'TTFSD','SAENZNOEE','DTFSD','OFF','HTFSD','OFF');
Fonnoee.ComputeTBDDS(); Fonnoee.ComputeTBDIS(); 
%AddStatFluctTBDIS(Fonnoee);
Fonnoee.DisplayTDBInfo();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Comparison With / Without FSD - Huge effect - 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot FSD on/off
figure(1)
hfsdoff = errorbar(Foff.qU-Foff.Q,Foff.TBDIS,Foff.TBDISE,...
    'sk','MarkerSize',4,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
errorbar_tick(hfsdoff,100);
PrettyFigureFormat
hold on
hfsdon = errorbar(Fon.qU-Fon.Q,Fon.TBDIS,Fon.TBDISE,...
    'ok','MarkerSize',4,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
errorbar_tick(hfsdon,100);
hold off
grid on
xlabel('qU-E_0 (eV)','FontSize',14);
ylabel('Counts','FontSize',14);
title(sprintf('KATRIN - %s - %g y - %g mcps',...
    Foff.TD,Foff.TimeYear,Foff.BKG_RateSec_i*1e3));
set(gca,'FontSize',12);
set(gca,'yscale','log');
    a = legend([hfsdoff hfsdon],...
        'FSD OFF','FSD ON (T-T, Saenz)','Location','NorthEast');
    legend(a,'boxoff');
publish_figure(1,'./figures/fsdonoff.eps')

% STEP 1: Ratio FSD On / Off
fsdratio  = Fon.TBDIS ./ Foff.TBDIS;
fsdratioe = (Fon.TBDISE.^2./Foff.TBDIS.^2 + Foff.TBDISE.^2./Fon.TBDIS.^2).^.5;
figure(2)
hfsdratio = errorbar(Foff.qU-Foff.Q,fsdratio,fsdratioe,...
    'sk','MarkerSize',6,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
errorbar_tick(hfsdratio,100);
PrettyFigureFormat
grid on
xlabel('qU-E_0 (eV)','FontSize',14);
ylabel('Ratio FSD ON/OFF','FontSize',14);
title(sprintf('KATRIN - %s - %g y - %g mcps',...
    Foff.TD,Foff.TimeYear,Foff.BKG_RateSec_i*1e3));
set(gca,'FontSize',12);
set(gca,'yscale','lin');
publish_figure(2,'./figures/fsdonoff-ratio.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Comparison With / Without FSD ONSET of Electronics excitations (EE)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot FSD with EE and without EE
figure(3)
hfsdon = errorbar(Fon.qU-Fon.Q,Fon.TBDIS,Fon.TBDISE,...
    'sk','MarkerSize',4,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
errorbar_tick(hfsdon,100);
PrettyFigureFormat
hold on
hfsdonnoee = errorbar(Fonnoee.qU-Fonnoee.Q,Fonnoee.TBDIS,Fonnoee.TBDISE,...
    'ok','MarkerSize',4,'MarkerFaceColor',.8*[1 1 1],'LineWidth',1);
errorbar_tick(hfsdonnoee,100);
hold off
grid on
xlabel('qU-E_0 (eV)','FontSize',14);
ylabel('Counts','FontSize',14);
title(sprintf('KATRIN - %s - %g y - %g mcps',...
    Fonnoee.TD,Fonnoee.TimeYear,Fonnoee.BKG_RateSec_i*1e3));
set(gca,'FontSize',12);
set(gca,'yscale','log');
    a = legend([hfsdon hfsdonnoee],...
        'FSD ON (T-T, Saenz)','FSD ON (T-T, Saenz) - No Electronics Excitations','Location','NorthEast');
    legend(a,'boxoff');
publish_figure(3,'./figures/fsdEEonoff.eps')

% STEP 2: Ratio FSD Electronics Excitations On / Off
fsdratio  = Fon.TBDIS ./ Fonnoee.TBDIS;
fsdratioe = (Fon.TBDISE.^2./Fonnoee.TBDIS.^2 + Fonnoee.TBDISE.^2./Fon.TBDIS.^2).^.5;
figure(4)
hfsdratio = errorbar(Fon.qU-Fon.Q,fsdratio,fsdratioe,...
    'sk','MarkerSize',6,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
errorbar_tick(hfsdratio,100);
PrettyFigureFormat
grid on
xlabel('qU-E_0 (eV)','FontSize',14);
ylabel('Ratio FSD EE ON/OFF - Integral','FontSize',14);
title(sprintf('KATRIN - %s - %g y - %g mcps',...
    Foff.TD,Foff.TimeYear,Foff.BKG_RateSec_i*1e3));
set(gca,'FontSize',12);
set(gca,'yscale','lin');
publish_figure(4,'./figures/fsdonoff-ratioint.eps')


%%% STEP 3: Ratio FSD Electronics Excitations On / Off
fsdratio  = Fon.TBDDS ./ Fonnoee.TBDDS;
fsdratioe = (Fon.TBDDS.^2./Fonnoee.TBDDS.^2 + Fonnoee.TBDDS.^2./Fon.TBDDS.^2).^.5;
figure(2)
hfsdratio = plot(Fon.Te-Fon.Q,fsdratio,...
    'sk','MarkerSize',2,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
PrettyFigureFormat
grid on
xlabel('T_e-E_0 (eV)','FontSize',14);
ylabel('Ratio FSD EE ON/OFF - Differential','FontSize',14);
title(sprintf('KATRIN - %s - %g y - %g mcps',...
    Foff.TD,Foff.TimeYear,Foff.BKG_RateSec_i*1e3));
set(gca,'FontSize',12);
set(gca,'yscale','lin');
publish_figure(4,'./figures/fsdonoff-ratiodiff.eps')

%% STEP 4: Changing the ES probability

% With FSD - With Electronics Excited States
clear Fon; 
Fon=InitKATRIN('TimeYear',TimeYear,'TTFSD','SAENZ','DTFSD','OFF','HTFSD','OFF');
Fon.ComputeTBDDS(); Fon.ComputeTBDIS(); 
%AddStatFluctTBDIS(Fon);
%Fon.DisplayTDBInfo();
fprintf(2,'Fon %g %g %g\n',Fon.TTNormGS,Fon.TTNormES,Fon.TTES_bias);

% With FSD - Without Electronics Excited States
clear Fonnoee;
Fonnoee=InitKATRIN('TimeYear',TimeYear,'TTFSD','SAENZNOEE','DTFSD','OFF','HTFSD','OFF');
Fonnoee.ComputeTBDDS(); Fonnoee.ComputeTBDIS(); 
%AddStatFluctTBDIS(Fonnoee);
%Fonnoee.DisplayTDBInfo();
fprintf(2,'Fonnoee %g %g %g\n',Fonnoee.TTNormGS,Fonnoee.TTNormES,Fonnoee.TTES_bias);

%%
figure(5)
fsdratio1 = plot(Fon.Te-Fon.Q,Fon.TBDDS ./ Fonnoee.TBDDS,'LineWidth',2,'Color','Blue');
hold on
C = linspecer(4);
for i=1:1:4
Fon.ComputeTBDDS('TTES_bias',-i/10); 
fprintf(2,'Fon %g %g %g\n',Fon.TTNormGS,Fon.TTNormES,Fon.TTES_bias);
hr(i) = plot(Fon.Te-Fon.Q,Fon.TBDDS./ Fonnoee.TBDDS,'LineWidth',2,'LineStyle','-','Color',C(i,:));
end
hold off
           
a = legend([fsdratio1 hr(1) hr(2) hr(3) hr(4)],'P_{EE}=0.42','P_{EE}=0.32','P_{EE}=0.22','P_{EE}=0.12','P_{EE}=0.0','Location','NorthEast') ; legend(a,'boxoff');
grid on
xlabel('T_e-E_0 (eV)','FontSize',14);
ylabel('Ratio FSD EE ON/OFF - Differential - Scale','FontSize',14);
title(sprintf('KATRIN - %s - %g y - %g mcps',...
    Foff.TD,Foff.TimeYear,Foff.BKG_RateSec_i*1e3),'FontSize',14);
set(gca,'FontSize',12);
set(gca,'yscale','lin');
PrettyFigureFormat
publish_figure(5,'./figures/fsdonoff-ratiodiff2.eps')

%%
figure(6)
Fon.ComputeTBDDS('TTES_bias',0); Fon.ComputeTBDIS(); 
fsdratio1 = plot(Fon.qU-Fon.Q,Fon.TBDIS ./ Fonnoee.TBDIS,'LineWidth',2,'Color','Blue');
hold on
C = linspecer(4);
for i=1:1:4
Fon.ComputeTBDDS('TTES_bias',-i/10); Fon.ComputeTBDIS(); 
fprintf(2,'Fon %g %g %g\n',Fon.TTNormGS,Fon.TTNormES,Fon.TTES_bias);
hr(i) = plot(Fon.qU-Fon.Q,Fon.TBDIS./ Fonnoee.TBDIS,'LineWidth',2,'LineStyle','-','Color',C(i,:));
end
hold off
           
a = legend([fsdratio1 hr(1) hr(2) hr(3) hr(4)],'P_{EE}=0.42','P_{EE}=0.32','P_{EE}=0.22','P_{EE}=0.12','P_{EE}=0.0','Location','NorthEast') ; legend(a,'boxoff');
grid on
xlabel('T_e-E_0 (eV)','FontSize',14);
ylabel('Ratio FSD EE ON/OFF - Integral - Scale','FontSize',14);
title(sprintf('KATRIN - %s - %g y - %g mcps',...
    Foff.TD,Foff.TimeYear,Foff.BKG_RateSec_i*1e3),'FontSize',14);
set(gca,'FontSize',12);
set(gca,'yscale','lin');
PrettyFigureFormat
publish_figure(6,'./figures/fsdonoff-ratiodintegral2.eps')

%%
figure(55)
    Fon.ComputeTBDDS('TTES_bias',0);
fsdratio1 = plot(Fon.Te-Fon.Q,Fon.TBDDS ./ Fonnoee.TBDDS,'LineWidth',2,'Color','Blue');
hold on
C = linspecer(4);
h1 = plot(Fon.TTexE_G,Fon.TTexP_G*100,'Color','Black','LineWidth',2);
hold on
h1e = plot(Fon.TTexE_E,Fon.TTexP_E*100,'Color','Blue','LineWidth',2,'LineStyle','-');
h2 = [0 0 0 0]; PE = [0 0 0 0];
for i=1:1:4
    Fon.ComputeTBDDS('TTES_bias',-i/10);
    fprintf(2,'Fon %g %g %g\n',Fon.TTNormGS,Fon.TTNormES,Fon.TTES_bias);
    h2(i) = plot(Fon.TTexE_E,Fon.TTexP_E*100,'Color',C(i,:),'LineWidth',2);
    PE(i) = Fon.TTNormES*100;
end
hold off
set(gca, 'YScale', 'log')  ;
set(gca, 'XScale', 'log');
grid on
axis([0.1 300 1e-2 50]);
xlabel('Excitation Energy','FontSize',12);
ylabel('Transition Probability (%)','FontSize',12);
title('T-T Ground and Excited States','FontSize',12)
set(gca,'FontSize',12);
strG = sprintf('Ground States: (set to P=%.2f %%)\n',Fon.TTNormGS*100);
strE = sprintf('Excited States: (set to P=%.2f %%)\n',100-Fon.TTNormGS*100);
strE1 = sprintf('Excited States: (set to P=%.2f %%)\n',PE(1));
strE2 = sprintf('Excited States: (set to P=%.2f %%)\n',PE(2));
strE3 = sprintf('Excited States: (set to P=%.2f %%)\n',PE(3));
strE4 = sprintf('Excited States: (set to P=%.2f %%)\n',PE(4));
a = legend([h1 h1e h2(1) h2(2) h2(3) h2(4)],strG,strE,strE1,strE2,strE3,strE4,'FontSize',10); 
legend(a,'boxoff','FontSize',10);
publish_figure(55,'./figures/TBD_FSD.eps')


%%
% Assign in base workspace
assignin('base', 'Fon', Fon);