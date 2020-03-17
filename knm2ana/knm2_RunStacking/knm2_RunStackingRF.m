% Run Stacking:closer look at response functions
% - single run response function vs. stacked run response function
%TwinOpt = {'Default','Twin_SameqUFlag','Twin_SameqUfracFlag','Twin_SameTime',...
%           'Twin_SameMTD','Twin_SameCDFlag','Twin_SameIsotopFlag','Sameall'};
TwinOpt = '';%'Twin_SameqUFlag';%'Twin_SameCDFlag';%'Twin_SameqUFlag';
if isempty(TwinOpt)
    TwinStr = '';
else
    TwinStr = ['_',TwinOpt];
end
RunListName = 'KNM2_Prompt';
TwinBias_Q = 18573.70;
savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_RunStackingRF_%s_E0%.2feV%s.mat',savedir,RunListName,TwinBias_Q,TwinStr);
%savename = sprintf('%stest.mat',savedir);

if exist(savename,'file')
    load(savename);
else
    RunAnaArg = {'RunList',RunListName,... % all KNM2 golden runs
        'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...
        'TwinBias_Q',TwinBias_Q,...
        'ROIFlag','14keV'};
    if ~isempty(TwinOpt)
        RunAnaArg = [RunAnaArg,TwinOpt,'ON'];
    end
    % build object of MultiRunAnalysis class
    A = MultiRunAnalysis(RunAnaArg{:});
    RFmodel = A.ModelObj.RF;
    
    RFsigma = std(A.SingleRunData.qU,0,2);
    A.ModelObj.MACE_Sigma = RFsigma;
    A.ModelObj.InitializeRF('RebinMode','Integral');
    RFmodelImp = A.ModelObj.RF;
    
%     A.ModelObj.MACE_Sigma = 2.*RFsigma; % enhance
%     A.ModelObj.InitializeRF('RebinMode','Integral');
%     RFmodelImp2 = A.ModelObj.RF;
%     
%     A.ModelObj.MACE_Sigma = 5.*RFsigma; % enhance
%     A.ModelObj.InitializeRF('RebinMode','Integral');
%     RFmodelImp5 = A.ModelObj.RF;
     
    A.ModelObj.MACE_Sigma = 10.*RFsigma; % enhance
    A.ModelObj.InitializeRF('RebinMode','Integral');
    RFmodelImp10 = A.ModelObj.RF;
    
    %reset
    A.ModelObj.MACE_Sigma = 0;
    A.ModelObj.InitializeRF;
    
    % load response functions of every model
    A.LoadSingleRunObj;
    TeModel = A.ModelObj.Te;
    qUModel = A.ModelObj.qU;
    
    nTe = 2250;
    Te = cell2mat(cellfun(@(x) x.Te(1:nTe),A.SingleRunObj,'UniformOutput',false)')';
    qU = cell2mat(cellfun(@(x) x.qU,A.SingleRunObj,'UniformOutput',false)')';
    
    RF = zeros(A.nRuns,nTe,A.ModelObj.nqU);
    for i=1:A.nRuns
        RF(i,:,:) = A.SingleRunObj{i}.RF(1:nTe,:);
    end
    
    nRuns = A.nRuns;
    RunList = A.RunList;
  
 %% interpolate RF, so that they all have same Te-qU
    RFinter= zeros(nRuns,numel(TeModel),numel(qUModel));
    for i=1:nRuns
        progressbar(i/nRuns);
        for q=1:numel(qUModel)
            RFinter(i,:,q) = squeeze(interp1(Te(i,:)-qU(i,q),RF(i,:,q),TeModel-qUModel(q),'linear','extrap'));
        end
    end
    
    RFinterMean = squeeze(mean(RFinter));
    RFinterStd = squeeze(std(RFinter));
   
    %% fit with with interpolated RF
    A.exclDataStart=A.GetexclDataStart(40);
    A.ModelObj.InitializeRF;
    A.fixPar = 'mNu E0 Bkg Norm';
    A.InitFitPar;
    A.Fit;
    FitResult_RFref = A.FitResult;
    
    A.ModelObj.RF = RFinterMean;
    A.Fit;
    FitResult_RFmean = A.FitResult;
    
    ISProb      = cell2mat(cellfun(@(x) x.ComputeISProb,A.SingleRunObj,'UniformOutput',false)');
    ISProbmodel = A.ModelObj.ComputeISProb;
    %% save
    save(savename,'Te','qU','RF','RFmodel','RFmodelImp','RFmodelImp10',...
                  'RFsigma','qUModel','TeModel','nRuns','RunList',...
                  'FitResult_RFmean','FitResult_RFref',...
                  'RFinter','RFinterMean','RFinterStd',...
                  'ISProb','RunAnaArg','ISProbmodel');
             
end

%% test start
A = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart=A.GetexclDataStart(40);
%new rf inter test
RFinter= zeros(nRuns,numel(TeModel),numel(qUModel));
for i=1:nRuns
    progressbar(i/nRuns);
    for q=1:numel(qUModel)
        RFinter(i,:,q) = squeeze(interp1(Te(i,:),RF(i,:,q),TeModel,'linear','extrap'));
    end
end
RFinterMean = squeeze(mean(RFinter));
RFinterStd = squeeze(std(RFinter));

A.Fit;
mNuSqref = A.FitResult.par(1);
A.ModelObj.RF = RFinterMean;
A.Fit;
mNuSqinter = A.FitResult.par(1);

% test end
%% correct for extrap mistakes
for i=1:numel(qUModel)
    RFinterMean(TeModel-qUModel(i)<=0,i) = 0;
end
%% plot 1: response function transmission edge area
pSingleObj = cell(nRuns,1);
qUi = 27;

f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
% for i=1:nRuns
%     pSingleObj{i} = plot(Te(i,:)-qU(i,qUi),RF(i,:,qUi),'LineWidth',1,'Color',rgb('PowderBlue'));
%     hold on;
% end
[l1,a1] = boundedline(TeModel-qUModel(qUi),RFinterMean(:,qUi),RFinterStd(:,qUi));
hold on;
l1.Color = rgb('DodgerBlue'); l1.LineWidth = 2;
a1.FaceColor = rgb('PowderBlue');
a1.FaceAlpha =0.5;
pref = plot(TeModel-qUModel(qUi),RFmodel(:,qUi),'--','Color',rgb('Orange'),'LineWidth',2);
pImp = plot(TeModel-qUModel(qUi),RFmodelImp(:,qUi),':','Color',rgb('Crimson'),'LineWidth',2);
%pImp2 = plot(TeModel-qUModel(qUi),RFmodelImp2(:,qUi),'-.','Color',rgb('DarkSlateGray'),'LineWidth',2);
%pImp5 = plot(TeModel-qUModel(qUi),RFmodelImp5(:,qUi),'-.','Color',rgb('Magenta'),'LineWidth',2);
pImp10 = plot(TeModel-qUModel(qUi),RFmodelImp10(:,qUi),'-.','Color',rgb('ForestGreen'),'LineWidth',2);
xlabel(sprintf('{\\itE} - %.0f (eV)',qUModel(qUi)));
ylabel('Transmission probability');
PrettyFigureFormat('FontSize',24)%pImp2,pImp5,
leg = legend([l1,a1,pref,pImp,pImp10],sprintf('Averaged RF     \\langleRF(x)\\rangle'),...
    sprintf('1\\sigma error band   \\langleRF(x)\\rangle'),...
    sprintf('Model average  RF(\\langlex\\rangle)'),...
    sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV   (\\sigma qU)',RFsigma(qUi)*1e3),...
   ... %sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV',2*RFsigma(qUi)*1e3),...
    ...%sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV',5*RFsigma(qUi)*1e3),...
    sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV (\\sigma qU \\times 10)',10*RFsigma(qUi)*1e3));
leg.EdgeColor = rgb('Silver');
leg.Location = 'southeast';
xlim([2.78,2.92])
saveplot1 = strrep(strrep(savename,'results','plots'),'.mat','_RFedge.pdf');
export_fig(f1,saveplot1);
fprintf('Save plot to %s \n',saveplot1);


%% fit with averaged (interpolated) response function
% result
fprintf('----------------------------------------\n')
fprintf('mNuSq = %.3f eV^2  (regular average RF) \n',FitResult_RFref.par(1))
fprintf('mNuSq =  %.3f eV^2  (averaged interp RF) \n',FitResult_RFmean.par(1))
fprintf('----------------------------------------\n')

%% plot2 : response function Oth scattering probabilities region
f2 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
pRF = plot(TeModel-qUModel(qUi),squeeze(RFinter(:,:,qUi)),':','LineWidth',1.5,'Color',...
             rgb('PowderBlue'));
hold on;
pmean = plot(TeModel-qUModel(qUi),RFinterMean(:,qUi),'-','Color',rgb('DodgerBlue'),'LineWidth',2);
pmodel = plot(TeModel-qUModel(qUi),RFmodel(:,qUi),'--','Color',rgb('Orange'),'LineWidth',2);
pISProbmin = plot(TeModel-qUModel(qUi),0.01*min(ISProb(1,:)).*ones(numel(TeModel),1),'-k','LineWidth',3);
pISProbmax = plot(TeModel-qUModel(qUi),0.01*max(ISProb(1,:)).*ones(numel(TeModel),1),'-k','LineWidth',3);
leg = legend([pRF(1),pmean,pmodel,pISProbmin],'Single runs response functions',...
              'Average response function',...
              'Stacked run response function',...
              'Single runs min/max 0th scattering probability',...
              'Location','northwest','EdgeColor',rgb('Silver'));
xlim([0,14]);
ylim([0.444,0.451])
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('{\\itE} - %.0f (eV)',qUModel(qUi)));
ylabel('Transmission probability');
saveplot2 = strrep(strrep(savename,'results','plots'),'.mat','_RF0Scat.pdf');
export_fig(f2,saveplot2);
fprintf('Save plot to %s \n',saveplot2);

%% plot 3 whole response functon
f3 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
[l1,a1] = boundedline(TeModel-qUModel(qUi),RFinterMean(:,qUi),50.*RFinterStd(:,qUi));
hold on;
l1.Color = rgb('DodgerBlue'); l1.LineWidth = 2.5;
a1.FaceColor = rgb('PowderBlue');
a1.FaceAlpha =0.5;
pref = plot(TeModel-qUModel(qUi),RFmodel(:,qUi),'--','Color',rgb('Orange'),'LineWidth',2);
xlabel(sprintf('{\\itE} - %.0f (eV)',qUModel(qUi)));
ylabel('Transmission probability');
PrettyFigureFormat('FontSize',24)%pImp2,pImp5,    
leg = legend([l1,a1,pref],sprintf('Average single runs RF  \\langleRF(x)\\rangle'),...
    sprintf('1\\sigma error band \\times 50        \\langleRF(x)\\rangle'),...
    sprintf('Stacked run RF              RF(\\langlex\\rangle)'),...
    sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV   (\\sigma qU)',RFsigma(qUi)*1e3),...
    sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV   (\\sigma qU \\times 10)',10*RFsigma(qUi)*1e3));
leg.EdgeColor = rgb('Silver');
leg.Location = 'southeast';
xlim([0,40])
ylim([0 1]);
 saveplot3 = strrep(strrep(savename,'results','plots'),'.mat','_RF.pdf');
 export_fig(f3,saveplot3);
 fprintf('Save plot to %s \n',saveplot3);
%% plot whole RF difference
% result
% def, def with const cross section, same CD: ~1e-3
% same qU: ~1e-7
% --> difference in response function comes from qU
f4 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
plot(TeModel-qUModel(qUi),zeros(numel(TeModel),1),'-k','LineWidth',2);
hold on;
y = RFinterMean(:,qUi)-RFmodel(:,qUi);
p1 = plot(TeModel-qUModel(qUi),y,'LineWidth',2,'Color',rgb('Orange'));
xlabel(sprintf('{\\itE} - %.0f (eV)',qUModel(qUi)));
ylabel('Transmission probability diff.');
PrettyFigureFormat('FontSize',24)%
leg = legend(p1,sprintf('Average single runs \\langleRF(x)\\rangle - Stacked run RF(\\langlex\\rangle)'));
leg.EdgeColor = rgb('Silver');
leg.Location = 'north';
xlim([-1,40])
%ylim([-6 7]*1e-4);
%ylim([min(y) max(y)]);
grid off

saveplot4 = strrep(strrep(savename,'results','plots'),'.mat','_RFdiff.pdf');
export_fig(f4,saveplot4);
fprintf('Save plot to %s \n',saveplot4);
%% scattering probabilities comparison at fixed E=18575
% result
% def, same qU, same qUfrac, def with const cross section:  have ~1e-5 differences
% same CD has no difference (1e-13 --> numerical) 
% --> influence of 1e-5 in scattering probabilities seems very small
% --> scattering probabilities are not the cause of difference in response function
f5 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
plot(-1:8,zeros(10,1),'-k','LineWidth',2);
hold on;
%plot(1:numel(ISProbmodel),ISProb-ISProbmodel,'o');
p1 = plot(0:numel(ISProbmodel)-1,mean(ISProb,2)-ISProbmodel,...
    ':o','LineWidth',2,'Color',rgb('Orange'),'MarkerFaceColor',rgb('Orange'),'MarkerSize',9);
PrettyFigureFormat('FontSize',24);
xlabel('Scattering');
ylabel(sprintf('Probability diff.'));
leg = legend(p1,sprintf('Average single runs - stacked run'));
leg.EdgeColor = rgb('Silver');
xlim([-0.3,7.3])
hold off;

saveplot5 = strrep(strrep(savename,'results','plots'),'.mat','_ScatteringProbs.pdf');
export_fig(f5,saveplot5);
fprintf('Save plot to %s \n',saveplot5);


