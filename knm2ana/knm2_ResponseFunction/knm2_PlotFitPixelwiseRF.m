% plot response function for KNM2 (data until 28th Oct,2019)
% % plot mean and std for pixelwise RF

% settings
RunList               = 'KNM2_Prompt';
exclDataStart         = 11; % 40eV range = 27 subruns
chi2                  = 'chi2Stat';
fixPar = 'E0 Bkg Norm'; % free Parameter !
pullFlag = 4;
% Init Model Object and covariance matrix object
Real = MultiRunAnalysis('RunList',RunList,'DataType','Real','AnaFlag','StackPixel','NonPoissonScaleFactor',1);

% Get pixelwise information
Real.ReadSingleRunData;Re
%% magnetic field in analyzing plane
Sigma_Ba_T = std(Real.SingleRunData.MACE_Ba_T);                                          % standard deviation over pixel for every run
Sigma_Ba_T_Rel = std(Real.SingleRunData.MACE_Ba_T)./mean(Real.SingleRunData.MACE_Ba_T);  % rel. standard deviation over pixel for every run
Mean_Ba_T =  mean(Real.SingleRunData.MACE_Ba_T'); % average Ba over pixels

meanSigma_Ba_T = mean(Sigma_Ba_T);            % mean standard deviation over pixels
meanSigma_Ba_T_Rel = mean(Sigma_Ba_T_Rel);    % mean rel.standard deviation over pixels

[plotHandle, cbHandle] = FPDViewer((Mean_Ba_T-mean(Mean_Ba_T))*1e6);
cbHandle.Label.String = sprintf('B_{ana} - \\langleB_{ana}\\rangle (10^{-6}T)');
cbHandle.FontSize = get(gca,'FontSize')+4;
cbHandle.Position(1) = 0.81;
savedirplot = strrep(savedir,'results','plots');
MakeDir(savedirplot);
savenameBa = [savedirplot,'knm2_PixelmapBa.pdf'];
export_fig(plotHandle,savenameBa,'-painters');
%% retarding poitential in analyzing plane: qUOffsets from runsummary
h5path = [getenv('SamakPath'), 'tritium-data/hdf5/',Real.DataSet,'/'];
h5name = ['RunSummary-Prompt3b-fpd00',num2str(Real.RunList(end)),'.h5'];
qUOffsets = h5read([h5path,h5name],'/RunSummary/RetardingEnergies/RetardingEnergyOffsets'); %V148
[plotHandle, cbHandle] = FPDViewer((qUOffsets-mean(qUOffsets))*1e3);
cbHandle.Label.String = sprintf('qU - \\langleqU\\rangle (mV)');
cbHandle.FontSize = get(gca,'FontSize')+4;
cbHandle.Position(1) = 0.81;
savenameqU = [savedirplot,'knm2_PixelmapqU.pdf'];
export_fig(plotHandle,savenameqU,'-painters');


%% retarding poitential in analyzing plane: qUOffsets from stacked data
% quick cross check, if pixel inhomogeniety stays the same when stacking the runs
% result: yes, it's the same
qU = permute(Real.SingleRunData.qU,[3,2,1]);
qURunStack = Real.StackWmean(Real.SingleRunData.qU,Real.SingleRunData.TimeperSubRunperPixel) ; %qU-values averaged over runs size(nqU,npixel)
subrun = 5;
meanqURunStack =mean(qURunStack,2);
[plotHandle, cbHandle] = FPDViewer((qURunStack(subrun,:)-meanqURunStack(subrun)).*1e3);
cbHandle.Label.String = sprintf('qU - \\langleqU\\rangle (mV)');
cbHandle.FontSize = get(gca,'FontSize')+4;
cbHandle.Position(1) = 0.81;

%% plot response function: mean and std for all active pixels
savedir = [getenv('SamakPath'),'knm2ana/knm2_ResponseFunction/results/'];
savename = [savedir,sprintf('knm2_ResponseFunction_Pixelwise2.mat')];
d = importdata(savename);

RFmean = mean(d.RF,3);  % average over active pixels
RFstd  = std(d.RF,0,3); % standard deviation over active pixels

subrun  = 15;
%[l1,a1] = boundedline(Real.ModelObj.Te-Real.ModelObj.qU(subrun),RFmean(:,subrun),RFstd(:,subrun));
%%
if ~isfield(d,'RFinterp')
 RFinterp = zeros(Real.ModelObj.nTe,numel(Real.PixList));
 for i=1:numel(Real.PixList)
    progressbar(i/numel(Real.PixList));
   %  plot(d.Te{i}-d.qU(subrun,i),squeeze(d.RF(:,subrun,i)))
   %  hold on;
   RFinterp(:,i) = interp1(d.Te{i}-d.qU(subrun,i),squeeze(d.RF(:,subrun,i)),Real.ModelObj.Te-Real.ModelObj.qU(subrun));
 end
 save(savename,'RFinterp','-append');
else
    RFinterp = d.RFinterp;
end
figRF = figure('Units','normalized','Position',[0.1,0.1,0.4,0.4]);
[l1,a1]  =  boundedline(Real.ModelObj.Te-Real.ModelObj.qU(subrun),mean(RFinterp,2),100.*std(RFinterp,0,2));
l1.LineWidth = 1; l1.Color = 'k';
a1.FaceColor = rgb('ForestGreen');
 PrettyFigureFormat('FontSize',24);
 xlabel('Surplus energy (eV)');
 ylabel('Transmission probability');
 leg = legend([l1,a1],'Average response function',sprintf('1\\sigma error band \\times 100')); legend boxoff; 
 leg.Location = 'northwest';
 xlim([-5 40]);
 savenameRF = [savedirplot,'knm2_PixelwiseRF.pdf'];
export_fig(figRF,savenameRF,'-painters');
%% use different RF and fit
Real = MultiRunAnalysis('RunList',RunList,'DataType','Real','AnaFlag','StackPixel','NonPoissonScaleFactor',1,'fixPar','mNu E0 Norm Bkg','exclDataStart',11); % re-init
TBDIS_i = Real.ModelObj.TBDIS;

Real.RunData.qU = d.qU_Uniform;
Real.RunData.TBDIS = d.TBDIS_Uniform;
Real.Fit;




