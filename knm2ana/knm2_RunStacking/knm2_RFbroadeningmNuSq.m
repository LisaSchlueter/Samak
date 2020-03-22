% ------------------------------------------------------------------------------
% Test broadening of response function
% - Neutrino mass bias as function of std of response/transmission function
% - Fit MC twins unsing different response function broadenings (RFsigma)
% - plot result
% - result: mNuSq bias follows approximately 2*sigma^2 formula
% January 20, Lisa
% ------------------------------------------------------------------------------

MACE_Sigma = 1e-3.*[0,1,5,10,15,20:10:100,125,150,200,250]; % standard deviation of gaussian broadening
% label
savedir  = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
savename = [savedir,sprintf('knm2_RFbroadeningmNuSq_MinSigma%.0fmeV_MaxSigma%.0fmeV.mat',...
               1e3*min(MACE_Sigma),1e3*max(MACE_Sigma))];

if exist(savename,'file')
    load(savename);
else
    
    % model settings
    range = 40;
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'DataType','Twin',...
        'fixPar','mNu E0 Bkg Norm',...
        'TwinBias_Q',18573.70,...
        'NonPoissonScaleFactor',1};
    
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    
    % init
    nSigma = numel(MACE_Sigma);
    mNuSq = zeros(nSigma,1);
    E0    = zeros(nSigma,1);
    RF    = zeros(nSigma,A.ModelObj.nTe,A.ModelObj.nqU);
    
    % reference model without broadening
    TBDIS = A.ModelObj.TBDIS; % Asimov spectra
    A.RunData.TBDIS = TBDIS;

    %% start broadening
    progressbar('RF broadening - calculate mNuSq bias');
    for i=1:nSigma
        progressbar(i/nSigma);
        A.ModelObj.MACE_Sigma = MACE_Sigma(i);
        A.ModelObj.InitializeRF('RebinMode','Integral');
        RF(i,:,:) = A.ModelObj.RF;
        A.Fit;
        mNuSq(i) = A.FitResult.par(1);
        E0(i)    = A.FitResult.par(2);
        
        % reset
        A.ModelObj.MACE_Sigma = 0;
        A.ModelObj.InitializeRF;
    end
    
    save(savename,'mNuSq','E0','RF','MACE_Sigma','range','TBDIS','RunAnaArg')
end

%% plot result

f123 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
p1 = plot(MACE_Sigma,-mNuSq,'LineWidth',3,'Color',rgb('DodgerBlue'));
hold on;
x = linspace(0,max(MACE_Sigma),100);
%p2 = plot(x,-(2*x.^2),'-.','LineWidth',3,'Color',rgb('Orange'));
xlabel(sprintf('\\sigma (eV)'));
ylabel(sprintf('\\Delta{\\itm}_\\nu^2 (eV^{ 2})'));
PrettyFigureFormat('FontSize',24); 
%leg = legend('Response function broadening',sprintf('\\Delta{\\itm}_\\nu^2 = -2\\sigma^2'));
%leg.Location = 'southwest';
%leg.EdgeColor = rgb('Silver');
ylim([min(-mNuSq),0.01])
title(sprintf('%.0f eV range',range),'FontWeight','normal','FontSize',get(gca,'FontSize'));


%save
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(f123,plotname);
fprintf('save plot to %s\n',plotname)


