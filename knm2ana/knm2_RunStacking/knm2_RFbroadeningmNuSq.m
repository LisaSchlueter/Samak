% ------------------------------------------------------------------------------
% Test broadening of response function
% - Neutrino mass bias as function of std of response/transmission function
% - Fit MC twins unsing different response function broadenings (RFsigma)
% - plot result
% - result: mNuSq bias follows approximately 2*sigma^2 formula
% January 20, Lisa
% ------------------------------------------------------------------------------

RFsigma = (0.0:0.05:0.7); % standard deviation of gaussian broadening

% label
savedir  = [getenv('SamakPath '),'knm2ana/knm2_RunStacking/results/'];
savename = [savedir,sprintf('knm2_RFbroadeningmNuSq_MinSigma%.0fmeV_MaxSigma%.0fmeV.mat',1e3*min(RFsigma),1e3*max(RFsigma))];

if exist(savename,'file')
    load(savename);
else
    RunList = 'KNM2_RW1';
    range = 40;
    RunAnaArg = {'RunList',RunList,...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Bkg Norm',...
        'TwinBias_Q',18573.70};
    %% set up model
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    %% reference model without broadening
    A.ModelObj.MACE_Sigma = 0;
    A.ModelObj.InitializeRF;
    A.Fit;
    mNuSq_ref = A.FitResult.par(1);
    E0_ref    = A.FitResult.par(2);
    
    % model with broadening
    mNuSq = zeros(numel(RFsigma),1);
    E0    = zeros(numel(RFsigma),1);
    
    progressbar('RF broadening - calculate mNuSq bias');
    for i=1:numel(RFsigma)
        progressbar(i/numel(RFsigma));
        A.ModelObj.MACE_Sigma = RFsigma(i);
        A.ModelObj.InitializeRF;
        A.Fit;
        A.ModelObj.MACE_Sigma = 0;
        A.ModelObj.InitializeRF;
        mNuSq(i) = A.FitResult.par(1);
        E0(i)    = A.FitResult.par(2);
    end
    
    save(savename,'mNuSq_ref','mNuSq','E0_ref','E0','RFsigma','range','RunAnaArg')
end

%% plot result

f123 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
p1 = plot(RFsigma*1e3,1e3.*(mNuSq-mNuSq_ref),'LineWidth',2.5,'Color',rgb('DodgerBlue'));
hold on;
p2 = plot(RFsigma*1e3,1e3.*(2*RFsigma.^2),'--','LineWidth',2.5,'Color',rgb('Orange'));
xlabel(sprintf('\\sigma_{RF} (meV)'));
ylabel(sprintf('\\Delta m_\\nu^2 (meV^2)'));
PrettyFigureFormat('FontSize',24);

leg = legend([p1,p2],'Response function broadening',sprintf('2 \\sigma_{RF}^2'));
leg.Location = 'northwest';
leg.EdgeColor = rgb('LightGray');
hold off;




