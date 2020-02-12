% estimate influence of T-minus decays on fake MC data (KNM2-like)
% Use real T-minus FSD
% Improved study with respect to knm2_TminusIons.m
% Lisa, November 2019

% settings
InitFile = @ref_FakeRun_KNM2_CD84_50days; % KNM2-like: 84% column density, 50 days, etc.
range = 40;
RecomputeFlag = 'OFF';

if range==90
    exclDataStart = 1;
elseif range == 40
    exclDataStart = 11;
end

WGTS_MolFrac_Tminus = linspace(1e-06,5*1e-03,10);%[1e-06,5e-06,1e-05,5e-05,1e-04,1e-04,1e-03,5e-03,1e-02];%(1e-06:1e-05:1e-04); % test fraction of Tminus ions
nTest = numel(WGTS_MolFrac_Tminus);

RunAnaArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FakeInitFile',InitFile,...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'exclDataStart',exclDataStart,... % 1==90 eV range ,11==40eV range (28 subruns)
    'chi2','chi2Stat',...
    'minuitOpt','min;migrad',...
    'NonPoissonScaleFactor',1,...
    'AnaFlag','StackPixel',...
    'fixPar','mNu E0 Norm Bkg'};

%labeling
savedir = [getenv('SamakPath'),'knm2ana/knm2_Plasma/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_Plasma_%.0eVrange_Tminus_Min%.2g_Max%.2g_nTest%.0f_FSDTminusSAENZ.mat',range,min(WGTS_MolFrac_Tminus),max(WGTS_MolFrac_Tminus),nTest)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename,'FitResult','TBDIS_i','TBDIS','RunAnaArg','qU','TimeSubrun','Bkg','WGTS_MolFrac_Tminus');
else
    %% model with T minus
    MT = RunAnalysis(RunAnaArg{:});
    MT.ModelObj.WGTS_MolFrac_Tm_i = 0; % reference spectrum without T minus
    MT.ModelObj.BKG_RateSec_i = 0.210; MT.ModelObj.SetFitBias(1);
    MT.ModelObj.TmFSD = 'SAENZ';
    MT.ModelObj.LoadFSD;
    MT.ModelObj.AdjustMolFrac;
    MT.ModelObj.ComputeTBDDS;
    MT.ModelObj.ComputeTBDIS;
    TBDIS_i = MT.ModelObj.TBDIS; % reference integral spectrum without T minus
    
    %% model without T minus
    M = RunAnalysis(RunAnaArg{:}); %model without T minus
    MT.ModelObj.TmFSD = 'SAENZ';
    MT.ModelObj.LoadFSD;
    MT.ModelObj.AdjustMolFrac;
    M.ModelObj.BKG_RateSec_i = 0.210; M.ModelObj.SetFitBias(1);
    %% compute and fit  test spectra with different t- fractions
    TBDIS = zeros(nTest,MT.ModelObj.nqU);
    FitResult    = cell(nTest,1);
    
    progressbar('T minus ions test');
    for i=1:nTest
        progressbar(i/nTest);
        MT.ModelObj.WGTS_MolFrac_Tm_i =  WGTS_MolFrac_Tminus(i); MT.ModelObj.SetFitBias(1);
        MT.ModelObj.AdjustMolFrac;
        MT.ModelObj.ComputeTBDDS;
        MT.ModelObj.ComputeTBDIS;
        TBDIS(i,:) = MT.ModelObj.TBDIS;
        
        M.RunData.TBDIS = TBDIS(i,:)';
        M.Fit;
        FitResult{i} = M.FitResult;
    end
    
    TimeSubrun = M.RunData.qUfrac.*M.RunData.TimeSec;
    qU = M.RunData.qU;
    Bkg = M.ModelObj.BKG_RateSec.*TimeSubrun;
    save(savename,'FitResult','TBDIS_i','TBDIS','RunAnaArg','qU','TimeSubrun','Bkg','WGTS_MolFrac_Tminus');
end
%% plots
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
%% plot 1: spectrum (sanity plot)
SanityPlot = 'ON';
if strcmp(SanityPlot,'ON')
    PlotIndex =nTest;
    
    WGTS_MolFrac_Tm = 0.02;
    % setup model with only T minus
    T =  RunAnalysis(RunAnaArg{:}); % Tminus model
    T.ModelObj.BKG_RateSec_i = 0.210; T.ModelObj.SetFitBias(1);
    TT_i = T.ModelObj.WGTS_MolFrac_TT;
    HT_i = T.ModelObj.WGTS_MolFrac_HT;
    DT_i = T.ModelObj.WGTS_MolFrac_DT;
    T.ModelObj.WGTS_MolFrac_DT = 0;%1e-15;%no DT
    T.ModelObj.WGTS_MolFrac_HT = 0;%1e-15;%no HT
    T.ModelObj.WGTS_MolFrac_TT = 0;%1e-15;%no TT
    
    T.ModelObj.TmFSD = 'SAENZ';
    T.ModelObj.LoadFSD;
    T.ModelObj.WGTS_MolFrac_Tm_i = WGTS_MolFrac_Tm; T.ModelObj.SetFitBias(1);%WGTS_MolFrac_Tminus(PlotIndex);
    T.ModelObj.AdjustMolFrac;
    T.ModelObj.ComputeTBDDS;
    T.ModelObj.ComputeTBDIS;
    TBDIS_T = T.ModelObj.TBDIS;
    
    figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
    pTT = plot(qU-18573.7,(TBDIS_i-Bkg)./TimeSubrun,'o-','LineWidth',2); % no T minus
    hold on;
    pTm = plot(qU-18573.7,TBDIS_T./TimeSubrun-T.ModelObj.BKG_RateSec,'o-','LineWidth',2); % only T minus
    PrettyFigureFormat('FontSize',24);
    leg = legend([pTT,pTm],...
        sprintf('%.3g%% T_2 , %.3g%% HT , %.3g%% of DT ',100*TT_i,100*HT_i,100*DT_i),...
        sprintf('%.3g%% T^{ -} only',100*WGTS_MolFrac_Tm));
    legend boxoff
    
    set(gca,'YScale','log');
    xlabel(sprintf('Retarding potential - 18573.7 (eV)'))
    ylabel('Rate (cps)');
    xlim([-90 10]);
    ylim([0.015 1e4]);
    
    % save
    savefile = [getenv('SamakPath'),sprintf('knm2ana/knm2_Plasma/plots/TminusSaenzFSD_spectrum_%.0fTm.pdf',WGTS_MolFrac_Tm)];
    export_fig(gcf,savefile,'-painters');
    fprintf('save plot to %s \n',savefile)
end

%% plot 2: result
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
Parameter = 'E0';%'mNu';
if strcmp(Parameter,'mNu')
    PlotPar = 1;
    ystr = sprintf('{\\it m}_\\nu^2 (eV^2)');
elseif strcmp(Parameter,'E0')
    PlotPar = 2;
    ystr = sprintf('{\\it E}_0^{fit} (meV)');
elseif strcmp(Parameter,'B')
    PlotPar = 3;
    ystr = sprintf('{\\it B (mcps)}');
elseif strcmp(Parameter,'N')
    PlotPar = 4;
    ystr = sprintf('{\\it N}');
end
y    = cell2mat(cellfun(@(x) x.par(PlotPar),FitResult,'UniformOutput',0));
yErr = cell2mat(cellfun(@(x) x.err(PlotPar),FitResult,'UniformOutput',0));

if ismember(Parameter,{'Nu','E0','B'})
    y = y.*1e3;
    yErr = yErr.*1e3;
end
x = WGTS_MolFrac_Tminus;
p = plot(x,y,'-o','LineWidth',2,'MarkerFaceColor',rgb('DodgerBlue'),'MarkerSize',8);
%p.CapSize = 0;
PrettyFigureFormat('FontSize',24);
xlabel('Fraction T^{ -}');
ylabel(ystr);
xlim([0,max(x).*1.05]);

saveplot2 = [plotdir,sprintf('knm2_TminosIonsFSD_%s.pdf',Parameter)];
export_fig(gcf,saveplot2,'-painters');
fprintf('save plot to %s \n',saveplot2)

