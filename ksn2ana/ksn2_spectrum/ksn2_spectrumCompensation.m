% illustration how positive m4 can be compensated

%% configure RunAnalysis object
%A = ref_KNM2_KATRIN_RegMTD;
R = RunAnalysis('RunNr',1,...
    'DataType','Fake',...
    'FakeInitFile',@ref_KNM2_KATRIN_RegMTD,...
    'fixPar','mNu E0 Norm Bkg');
R.exclDataStart = R.GetexclDataStart(40);
M = R.ModelObj;
M.BKG_RateSec_i =0.22;
M.normFit_i     = 1;
M.ComputeTBDDS;
M.ComputeTBDIS;
R.RunData.TBDIS = M.TBDIS;

Time = M.TimeSec.*M.qUfrac;
Data = M.TBDIS./Time;


%% sterile + active nu-mass = 0
mNu4Sq = 50;
sin2T4  = 0.4;
M.SetFitBiasSterile(mNu4Sq,sin2T4)
M.ComputeTBDDS;
M.ComputeTBDIS;
Model_Sterile = M.TBDIS./Time; % 60% active with m_nu = 0, 40% sterile

%% plot 1
GetFigure;

pCombi = plot(M.qU-M.Q_i,Model_Sterile./Data,'-','LineWidth',2);
hold on;
pSterile = plot(M.qU-M.Q_i,(Model_Sterile-0.6.*Data)./Data,':','LineWidth',2);
pActive  = plot(M.qU-M.Q_i,(0.6.*Data)./Data,'-.','LineWidth',2);
pd = plot(M.qU-M.Q_i,Data./Data,'.k','MarkerSize',12);
PrettyFigureFormat;
xlabel(sprintf('Retarding potential - {\\itE}_0 (eV)'));
ylabel('Ratio');
xlim([-42 42]);
ylim([0.3 1.1])
leg = legend([pd,pActive,pSterile,pCombi],'MC data',...
    sprintf('Active branch {\\itm}_\\nu^2 = 0 eV^2'),...
    sprintf('Sterile branch {\\itm}_4^2= %.2g eV^2 , |{\\itU}_{e4}|^2 = %.3g',mNu4Sq,sin2T4),...
    'Active + Sterile branch','Location','east');
PrettyLegendFormat(leg);
%%
R.Fit;
mNuSqFit = R.FitResult.par(1);

%Model_SterilePlusNuMass = M.TBDIS./Time;

%%
M.SetFitBias(0);
M.SetFitBiasSterile(mNu4Sq,sin2T4)
M.mnuSq_i = mNuSqFit;
M.ComputeTBDDS;
M.ComputeTBDIS;
Model_NuMass = M.TBDIS./Time;

%% plot 2
GetFigure;
pCombi = plot(M.qU-M.Q_i,Model_NuMass./Data,'-','LineWidth',2);
hold on;
pSterile = plot(M.qU-M.Q_i,(Model_Sterile-0.6.*Data)./Data,':','LineWidth',2);
pActive  = plot(M.qU-M.Q_i,(Model_NuMass-(Model_Sterile-0.6.*Data))./Data,'-.','LineWidth',2);
pd = plot(M.qU-M.Q_i,Data./Data,'.k','MarkerSize',12);
PrettyFigureFormat;
xlabel(sprintf('Retarding potential - {\\itE}_0 (eV)'));
ylabel('Ratio');
xlim([-42 42]);
ylim([0.3 1.1])
leg = legend([pd,pActive,pSterile,pCombi],'MC data',...
    sprintf('Active branch {\\itm}_\\nu^2 = %.3g eV^2',mNuSqFit),...
    sprintf('Sterile {\\itm}_4^2= %.2g eV^2 , |{\\itU}_{e4}|^2 = %.3g',mNu4Sq,sin2T4),...
    'Active + Sterile branch','Location','east');
PrettyLegendFormat(leg);


