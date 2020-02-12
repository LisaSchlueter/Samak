% Test fit to Tminus fraction on fake MC data


% settings
InitFile = @ref_FakeRun_KNM2_CD84_50days; % KNM2-like: 84% column density, 50 days, etc.
range = 90;
RecomputeFlag = 'ON';

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

%% create MC data with Tminus fraction
D = RunAnalysis(RunAnaArg{:}); 
D.ModelObj.WGTS_MolFrac_Tm_i = 2*1e-03;
D.ModelObj.BKG_RateSec_i = 0.210; D.ModelObj.SetFitBias(1);
D.ModelObj.TmFSD = 'SAENZ';
D.ModelObj.LoadFSD;
D.ModelObj.AdjustMolFrac;
D.ModelObj.ComputeTBDDS;
D.ModelObj.ComputeTBDIS;
TBDIS = D.ModelObj.TBDIS; % MC spectrum with T minus

%% model without T minus
M = RunAnalysis(RunAnaArg{:});
M.ModelObj.BKG_RateSec_i = 0.210; M.ModelObj.SetFitBias(1);
M.ModelObj.TmFSD = 'SAENZ';
M.ModelObj.LoadFSD;

%% Fit MC data
M.RunData.TBDIS = TBDIS;
M.fixPar = 'mNu E0 Norm Bkg FracTm'; M.InitFitPar; % additional free fit parameter: fractio of  T- 
M.Fit;




