function Fit_numass(varargin)
%
%            FIT INTEGRAL SPECTRUM
%           Fit Active Neutrino Mass
%
%          Th. Lasserre - CEA Saclay
%                February 2017
%

clear par ; clear parnomix;
addpath(genpath('../../../Samak2.0'));

Nfit    = 1; % Number of Fits

% Parametrization: True Value
mnuSq_t = (0.5)^2;

% Loop on fits
npar       = 6;
fit_p      = ones(npar,Nfit);
fit_xi2min = ones(1,Nfit);

% Parser
p = inputParser;

p.addParameter('nTeBinningFactor',50,@(x)isfloat(x) && x>0);

% KATRIN GENERAL SETTINGS
p.addParameter('TimeSec',3*86500*365,@(x)isfloat(x) && x>0);
p.addParameter('WGTS_CD_MolPerCm2',5e17);
p.addParameter('WGTS_B_T',2.52);
p.addParameter('MACE_Bmax_T',4.2);
p.addParameter('MACE_Ba_T',6e-4);

p.addParameter('Mode','Read',@(x)ismember(x,{'DataTBD', 'Sim', 'Read'}));
p.addParameter('TD','Flat30');
p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
p.addParameter('BKG_RateAllFPDSec',0.3,@(x)isfloat(x) && x>=0);
p.addParameter('mnuSq_i',0,@(x)isfloat(x));

% WGTS: Flags for FSD: T-T / D-T / H-T
p.addParameter('DTFSD','DOSS',@(x)ismember(x,{'OFF','DOSS'}));
p.addParameter('T2purity',0.99,@(x)isfloat(x) && x>0);

% Flag for Doppler Effect 
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','numConv','matConv'}));

% Transmission Function
p.addParameter('KTFFlag','Compute',@(x)ismember(x,{'OFF','SSC',...
    'MACE','MACE+WGTSIS','TL','SG','Kle',...
    'SSCW_DCperqU','SSC_BC','Compute'}));
p.addParameter('recomputeRF','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

% Parameters: Calculation precision
nTeBinningFactor = p.Results.nTeBinningFactor;

% Parameters: KATRIN GENERAL SETTINGS
Mode              = p.Results.Mode;
TD                = p.Results.TD;
WGTS_CD_MolPerCm2 = p.Results.WGTS_CD_MolPerCm2;
WGTS_B_T          = p.Results.WGTS_B_T;
MACE_Bmax_T       = p.Results.MACE_Bmax_T;
MACE_Ba_T         = p.Results.MACE_Ba_T; 

BKG_Flag          = p.Results.BKG_Flag;
BKG_Type          = p.Results.BKG_Type;
BKG_RateAllFPDSec = p.Results.BKG_RateAllFPDSec;
mnuSq_i           = p.Results.mnuSq_i;
TimeSec           = p.Results.TimeSec;

% Doppler Effect
DopplerEffectFlag = p.Results.DopplerEffectFlag;

% TBD: Flag FSD's
DTFSD           = p.Results.DTFSD;
T2purity        = p.Results.T2purity;

% Transmission Function 
KTFFlag         = p.Results.KTFFlag;
recomputeRF     = p.Results.recomputeRF;

% Default
opt_calc = {...
    'nTeBinningFactor',nTeBinningFactor};

opt_katrin = {...
    'Mode',Mode,...
    'TD',TD,...
    'TimeSec',TimeSec,...
    'mnuSq_i',mnuSq_i};

opt_wgts = {...
    'WGTS_CosMaxAAngle',0.630570414481244,...
    'WGTS_Tp',T2purity,...
    'WGTS_MolFrac_TT',0.95,...
    'WGTS_MolFrac_DT',0,...
    'WGTS_MolFrac_HT',0,...
    'WGTS_MolFracRelErr_TT',0,...
    'WGTS_MolFracRelErr_DT',0,...
    'WGTS_MolFracRelErr_HT',0,...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2,...
    'WGTS_B_T',WGTS_B_T,...
    'WGTS_Temp',30};

opt_mace = {...
    'MACE_Bmax_T',MACE_Bmax_T,...
    'MACE_Ba_T',MACE_Ba_T};

opt_wgtsmace = {...
    'KTFFlag','SSC_BC',...
    'recomputeRF',recomputeRF};

opt_fpd = {...
    'FPD_Segmentation','OFF','FPD_Pixel',2};

opt_bkg = {...
    'BKG_Flag',BKG_Flag,...
    'BKG_Type',BKG_Type,...
    'BKG_RateAllFPDSec',BKG_RateAllFPDSec};

opt_fsd= {...
    'TTFSD','DOSS',...
    'DTFSD','OFF',...
    'HTFSD','OFF'};

opt_doppler = {'DopplerEffectFlag',DopplerEffectFlag};

opt_integration = {'IStype','SIMPFAST'};

% Tritium spectrum definition
A = TBD(...
    opt_calc{:},...
    opt_katrin{:},...
    opt_wgts{:},...
    opt_mace{:},...
    opt_fpd{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_bkg{:},...
    opt_doppler{:},...
    opt_integration{:},...
    'mnuSq_i',mnuSq_t);

%disp(A.qU');

% Init
epsilon = 1e-99;
i_mnu = 0;
i_Q = 0;
i_B = 0;
i_N = 0;
i_m4 = 0;
i_s4 = 0;

% Number of Fits
tic
for f=1:1:Nfit
    
    % Initializing Fit
    A.ComputeTBDDS(); A.ComputeTBDIS(); AddStatFluctTBDIS(A);
    Data = [A.qU,A.TBDIS,A.TBDISE];
    DataTBD = {Data,A};
    parnames = ['mSq dQ B dN m4 s4'];
    ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4];
    tmparg = sprintf(['set pri -1;'...
        'fix  5 6 ;'...
        'set now; minos ; imp' ]);
    Args = {ParIni, DataTBD, '-c',tmparg};
    [par, err, chi2min, errmat] = fminuit('NuMassChi2',Args{:});
    fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
    
    fprintf('===============================================\n');
    fprintf('  m^2       = %g ± %g\n', (A.mnuSq_i+par(1)),err(1));
    fprintf('  dQ        = %g ± %g\n', (par(2)),err(2));
    fprintf('  B         = %g ± %g\n', (A.BKG_RateSec_i+par(3)),err(3));
    fprintf('  dN        = %g ± %g\n', par(4),err(4));
    fprintf('  Chi2/dof  = %g/%g\n', chi2min,A.nqU);
    fprintf('==========================================\n');
    
end
toc 

%% Plot Results
fig1 = figure(1)
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,1]);
subplot(2,1,1)
hfit = plot(Data(:,1)-A.Q_i,NuMassModel4par(par,A)./(A.qUfrac.*A.TimeSec),...
'LineWidth',1,'LineStyle','-','Color','Black');
hold on;
hdata = errorbar(Data(:,1)-A.Q_i,Data(:,2)./(A.qUfrac.*A.TimeSec),Data(:,3)./(A.qUfrac.*A.TimeSec),...
     'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
%errorbar_tick(hdata,200);
hold off
grid on
xlabel('qU-E_0 (eV)','FontSize',14);
ylabel('Rate (cps)','FontSize',14);
%title(sprintf('KATRIN - %s - %g y - %g mcps',A.TD,A.TimeSec/86400/365,A.BKG_RateSec_i));
set(gca,'FontSize',12);
set(gca,'yscale','log');
mydata = sprintf('Simulation: m_{eff}=%.2f eV \n',sqrt(abs(A.mnuSq_i)));
myfit = sprintf('Fit: m_{eff}=%.2f \\pm %.2f eV',sqrt(abs(A.mnuSq_i+par(1))),err(1));
mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-3);
a = legend([hdata hfit ],mydata,myfit,mychi2,'Location','NorthEast') ; 
axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min(Data(:,2)./(A.qUfrac.*A.TimeSec)) max(Data(:,2)./(A.qUfrac.*A.TimeSec))*1.1])
PrettyFigureFormat
 
subplot(2,1,2)
parnomix = zeros(1,4,1); parnomix = par; parnomix(1) = -A.mnuSq;
hfit = plot(Data(:,1)-A.Q_i,...
    NuMassModel4par(par,A)./NuMassModel4par(parnomix,A),...
    'Color','Black','LineWidth',1,'LineStyle','-');
hold on;
hdata = errorbar(Data(:,1)-A.Q,...
      Data(:,2)./NuMassModel4par(parnomix,A),...
      Data(:,3)./NuMassModel4par(parnomix,A),...
     'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
%errorbar_tick(hdata,200);
hold off;
grid on
xlabel('qU-E_0 (eV)','FontSize',14);
ylabel('Simulation/Zero-mass','FontSize',14);
set(gca,'FontSize',12);
a= legend([hdata hfit],'Simulation 3y - B=300 mcps','Fit','Location','NorthWest');
%legend(a,'boxoff');
%axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min((Data(:,2)./NuMassModel4par(parnomix,A))) max(Data(:,3)./NuMassModel4par(parnomix,A))*3])
xlim([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1]);
myname = sprintf('fitoutput'); publish_figure(1,myname);
PrettyFigureFormat
end