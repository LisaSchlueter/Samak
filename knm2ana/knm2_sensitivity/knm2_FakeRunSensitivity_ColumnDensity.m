% Compute neutrino mass sensitivity before data taking with fake MC run
% different measurement times
% sensitivity as a function of column density
% 1st Oct 19, Lisa
TimeDays = 50;
savedir = [getenv('SamakPath'),'knm2ana/knm2_sensitivity/results/'];
MakeDir(savedir);

switch TimeDays
    case 60
        savefile = [savedir,'knm2_FakeRunSensitivity_ColumnDensity_60days.mat'];
        InitFile = @ref_FakeRun_KNM2_CD100; %100 column density, 60 days
    case 1000
        savefile = [savedir,'knm2_FakeRunSensitivity_ColumnDensity_1000days.mat'];
        InitFile = @ref_FakeRun_KNM2_CD100_1000days; %100 column density, 1000 days
    case 50
        savefile = [savedir,'knm2_FakeRunSensitivity_ColumnDensity_50days.mat'];
        InitFile = @ref_FakeRun_KNM2_CD100_50days; %100 column density, 1000 days
end

if ~exist(savefile,'file')
CommonArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','Sibille0p5eV',...
    'fixPar','5 6 7 8 9 10 11',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'exclDataStart',12,... % 40eV range (28 subruns)
    'chi2','chi2Stat'};
R = RunAnalysis(CommonArg{:},'FakeInitFile',InitFile);
R.Fit;
B_bias = R.FitResult.par(3);
N_bias = R.FitResult.par(4);
%%
WGTS_CD = 5*1e17.*[(0.1:0.1:0.9),0.95,1];                 % column density
mNuSqErr = zeros(numel(WGTS_CD),1);                       % sensitivity
NormTBDDS = zeros(numel(WGTS_CD),1);                      % differential spectrum norm factor (sanity check)
NumElectrons =zeros(numel(WGTS_CD),1);                    % number of electrons in fit range (sanity check)
RF = zeros(numel(WGTS_CD),R.ModelObj.nTe,R.ModelObj.nqU); % response function (sanity check)

for i=1:numel(WGTS_CD)
    R.ModelObj.WGTS_CD_MolPerCm2 = WGTS_CD(i);
    R.ModelObj.AdjustRF;
    
    R.ModelObj.ComputeTBDDS('B_bias',B_bias,'N_bias',N_bias); R.ModelObj.ComputeTBDIS;
    R.RunData.TBDIS  = R.ModelObj.TBDIS;
    R.RunData.TBDISE = R.ModelObj.TBDISE;
    
    R.Fit;
    
    mNuSqErr(i) = R.FitResult.err(1);
    NormTBDDS(i) = R.ModelObj.NormFactorTBDDS;
    NumElectrons(i) = sum(R.ModelObj.TBDIS(R.exclDataStart:end));
    RF(i,:,:) = R.ModelObj.RF;
    
    R.ModelObj.SetFitBias(0); %reset fit init
end

save(savefile,'mNuSqErr','WGTS_CD','CommonArg','InitFile','NormTBDDS','NumElectrons','RF');

%% add 84% (KNM2 setting)
R.ModelObj.WGTS_CD_MolPerCm2 =0.84*5*1e17;
R.ModelObj.AdjustRF;

R.ModelObj.ComputeTBDDS('B_bias',B_bias,'N_bias',N_bias); R.ModelObj.ComputeTBDIS;
R.RunData.TBDIS  = R.ModelObj.TBDIS;
R.RunData.TBDISE = R.ModelObj.TBDISE;

R.Fit;
mNuSqErrKNM2 = R.FitResult.err(1); 
R.ModelObj.SetFitBias(0); %reset fit init

save(savefile,'mNuSqErrKNM2','-append');
else
    load(savefile);
end

%% sanity checks: number of electrons, normalization factor, ...
fprintf('---------------------------------------------------- \n');
fprintf('Sanity check \n');
fprintf('%.0f%% column density: NormFactorTBDDS= %.3g , Electrons = %.0f \n',...
    [WGTS_CD/(5*1e17)*100';NormTBDDS';NumElectrons']);
fprintf('---------------------------------------------------- \n');
fprintf('Results \n');
fprintf('%.0f%% column density: sigma(m_nu^2) = %.3g eV^2 \n',...
    [WGTS_CD/(5*1e17)*100';mNuSqErr']);
fprintf(2,'%.0f%% column density: sigma(m_nu^2) = %.3g eV^2 \n',...
    [84';mNuSqErrKNM2']);
fprintf('---------------------------------------------------- \n');
%% display
fig1 =figure(1);
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.4, 0.4]);
PlotIndex = mNuSqErr>0.01;
plotx = linspace(min(WGTS_CD(PlotIndex)),max(WGTS_CD(PlotIndex)),100)./(5*1e17)*100;
ploty = interp1(WGTS_CD(PlotIndex)./(5*1e17)*100,mNuSqErr(PlotIndex),plotx,'spline');

p1 =plot(plotx,ploty,...%WGTS_CD(PlotIndex)./(5*1e17)*100,mNuSqErr(PlotIndex),...
    '-','LineWidth',4,'Color',rgb('DarkMagenta'));
hold on;
e1 =errorbar(WGTS_CD(PlotIndex)./(5*1e17)*100,mNuSqErr(PlotIndex),zeros(sum(PlotIndex),1),...
    'o','LineWidth',4,'MarkerFaceColor',rgb('DarkMagenta'),'MarkerEdgeColor',rgb('DarkMagenta'));
e1.CapSize = 0;
PrettyFigureFormat('FontSize',24);
ylim([0.95*min(mNuSqErr),max(mNuSqErr)]);
xlim([min(WGTS_CD/(5*1e17)*100)-2,max(WGTS_CD/(5*1e17)*100)+2]);
p2 =plot(84.*[1,1],[min(ylim),max(ylim)],'-','LineWidth',p1.LineWidth,'Color',rgb('Silver'));
xlabel('Column density (%)');
ylabel(sprintf('1\\sigma sensitivity on m_\\nu^2'));
leg = legend([p1,p2],sprintf('Sensitivity %.0f days (KNM2 settings)',TimeDays),...
    [sprintf('\\sigma(m_\\nu^2) = %.2f eV',mNuSqErrKNM2),'^{ 2} @ 84% column density']);
legend boxoff;
leg.Location='north';
export_fig(fig1,strrep(strrep(savefile,'results','plots'),'.mat','.pdf'),'-painters');
hold off;


