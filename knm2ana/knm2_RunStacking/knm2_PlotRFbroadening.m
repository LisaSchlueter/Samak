% plot response functions with and without broadening

% settings
RecomputeFlag = 'OFF';
RFsigma = 0.2;       % standard deviation of gaussian broadening
RebinMode = 'Fast';  % Fast (good for RFsigma >= 0.1eV) or Integral (more precise, but slow)

savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
savename = [savedir,sprintf('knm2_BroadenRF_sigma%.0fmeV_%s.mat',1e3*RFsigma,RebinMode)];
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    RunList = 'KNM2_RW1';
    RunAnaArg = {'RunList',RunList,...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Bkg Norm',...
        'TwinBias_Q',18573.70};
    A = MultiRunAnalysis(RunAnaArg{:});
    
    % reference model without broadening
    A.ModelObj.MACE_Sigma = 0;
    A.ModelObj.InitializeRF;
    RF_ref = A.ModelObj.RF;
    
    % model with broadening
    A.ModelObj.MACE_Sigma = RFsigma;
    A.ModelObj.InitializeRF('RebinMode',RebinMode);
    RF_broadened = A.ModelObj.RF;
    
    qU = A.ModelObj.qU;
    Te = A.ModelObj.Te;
    
    save(savename,'qU','Te','RF_ref','RF_broadened','RunAnaArg');
end

%% plot response functions (regular + broadened)
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
qUindex = 20; % retaring potential (index)
pref = plot(Te-qU(qUindex),RF_ref(:,qUindex),'LineWidth',2.5);
hold on;
pb = plot(Te-qU(qUindex),RF_broadened(:,qUindex),'-','LineWidth',pref.LineWidth);
leg = legend(sprintf('\\sigma_{RF} = 0 meV'),sprintf('\\sigma_{RF} = %.0f meV',RFsigma*1e3),'Location','northwest');
leg.EdgeColor = rgb('LightGray');
xlabel('Surplus energy (eV)');
ylabel('Transmission probability');
PrettyFigureFormat('FontSize',24);
xlim([-1,6]);
ylim([-0.02 0.47]);
hold off;


