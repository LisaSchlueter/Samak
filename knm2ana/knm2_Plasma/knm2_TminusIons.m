% MC study
% Calculate effect of T-minus ions on neutrino mass
% no FSDs considered
% T-minus as T_2 without FSD but shitef endpoint (very approximate)
% more accurate study will follow
% Lisa, November 2019

% settings
InitFile = @ref_FakeRun_KNM2_CD84_50days; % KNM2-like: 84% column density, 50 days, etc.
range = 40;
RecomputeFlag = 'ON';

if range==90
    exclDataStart = 1;
elseif range == 40
    exclDataStart = 11;
end

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

WGTS_MolFrac_Tminus = linspace(1e-06,2*1e-03,10);%[1e-06,5e-06,1e-05,5e-05,1e-04,1e-04,1e-03,5e-03,1e-02];%(1e-06:1e-05:1e-04); % test fraction of Tminus ions
nTest = numel(WGTS_MolFrac_Tminus);
FSDFlag = 'OFF';%BlindingKNM2'; %BlindingKNM2 %FSD used for Tminus ion

%labeling
savedir = [getenv('SamakPath'),'knm2ana/knm2_Plasma/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_Plasma_%.0eVrange_Tminus_Min%.2g_Max%.2g_nTest%.0f_FSD-%s.mat',range,min(WGTS_MolFrac_Tminus),max(WGTS_MolFrac_Tminus),nTest,FSDFlag)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename,'FitResult','TBDIS_Tminus','TBDIS_TT','TimeSubrun','qU','Bkg','WGTS_MolFrac_TT_i');
else
    %% set up model, compute fake run if necessary
    M_TT = RunAnalysis(RunAnaArg{:}); %normal TT
    WGTS_MolFrac_TT_i = M_TT.ModelObj.WGTS_MolFrac_TT;
    
    M_Tm =  RunAnalysis(RunAnaArg{:}); % Tminus model
    M_Tm.ModelObj.TTFSD = FSDFlag; %no FSD -> just different endpoint
    M_Tm.ModelObj.DTFSD = FSDFlag;
    M_Tm.ModelObj.HTFSD = FSDFlag;
    M_Tm.ModelObj.WGTS_MolFrac_DT = 1e-15;%no DT
    M_Tm.ModelObj.WGTS_MolFrac_HT = 1e-15;%no HT
    M_Tm.ModelObj.Q_i = M_TT.ModelObj.Q_i + 7.45; % to do: what is the endpoint?
    M_Tm.ModelObj.AdjustMolFrac;
    
    %% Compute integral spectra with different fractions of Tminus
    TBDIS_TT     = zeros(nTest,M_TT.ModelObj.nqU);
    TBDIS_Tminus = zeros(nTest,M_TT.ModelObj.nqU);
    FitResult    = cell(nTest,1);
    
    for i=1:nTest
        % Model with reduced TT
        M_TT.ModelObj.WGTS_MolFrac_TT = WGTS_MolFrac_TT_i-WGTS_MolFrac_Tminus(i);
        M_TT.ModelObj.AdjustMolFrac; % Adjust normalization factor and FSD
        M_TT.ModelObj.ComputeTBDDS;
        M_TT.ModelObj.ComputeTBDIS;
        TBDIS_TT(i,:) = M_TT.ModelObj.TBDIS;
        
        % Model with small amout of Tminus
        M_Tm.ModelObj.WGTS_MolFrac_TT = WGTS_MolFrac_Tminus(i);
        M_Tm.ModelObj.AdjustMolFrac; % Adjust normalization factor and FSD
        M_Tm.ModelObj.ComputeTBDDS;
        M_Tm.ModelObj.ComputeTBDIS;
        TBDIS_Tminus(i,:) = M_Tm.ModelObj.TBDIS;
    end
    
    %% Fit spectra
    % re-init model assuming no Tminus
    M = RunAnalysis(RunAnaArg{:}); %normal TT
    
    TBDIS = TBDIS_TT + TBDIS_Tminus; % sum: TT and Tminus
    
    progressbar('T minus ions test');
    for i=1:nTest
        progressbar(i/nTest);
        M.RunData.TBDIS = TBDIS(i,:)';
        M.Fit;
        FitResult{i} = M.FitResult;
    end
    
    TimeSubrun = M_TT.RunData.qUfrac.*M_TT.RunData.TimeSec;
    qU = M_TT.RunData.qU;
    Bkg = M_TT.ModelObj.BKG_RateSec.*TimeSubrun;
    save(savename,'FitResult','TBDIS_Tminus','TBDIS_TT','RunAnaArg','qU','TimeSubrun','Bkg','WGTS_MolFrac_TT_i');
end

%% plots
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);

%% Sanity plot: spectra
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
PlotIndex =nTest;

pTT = plot(qU-18573.7,(TBDIS_TT(PlotIndex,:)'-Bkg)./TimeSubrun,'o-','LineWidth',2);
hold on;
pTm = plot(qU-18573.7,(TBDIS_Tminus(PlotIndex,:)'-Bkg)./TimeSubrun,'o-','LineWidth',pTT.LineWidth);
leg = legend([pTT,pTm],...
    sprintf('%.3g%% T_2',100*(WGTS_MolFrac_TT_i-WGTS_MolFrac_Tminus(PlotIndex))),...
    sprintf('%.3g%% T^{ -}',100*WGTS_MolFrac_Tminus(PlotIndex)));
legend boxoff
PrettyFigureFormat('FontSize',24);
leg.FontSize = get(gca,'FontSize');
set(gca,'YScale','log');
xlabel(sprintf('Retarding potential - 18573.7 (eV)'))
ylabel('Rate (cps)');
xlim([-90 10]);
ylim([9.9e-06 1e5]);
saveplot1 = [plotdir,sprintf('knm2_TminosIons_Spectra_%.2gTm.pdf',WGTS_MolFrac_Tminus(PlotIndex))];
export_fig(gcf,saveplot1,'-painters');
%% plot
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
Parameter = 'mNu';%'mNu';
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
xlabel('Molecular fraction T^{ -}');
ylabel(ystr);
xlim([0,max(x).*1.05]);

saveplot2 = [plotdir,sprintf('knm2_TminosIons_%s_%.2gTm.pdf',Parameter,WGTS_MolFrac_Tminus(PlotIndex))];
export_fig(gcf,saveplot2,'-painters');