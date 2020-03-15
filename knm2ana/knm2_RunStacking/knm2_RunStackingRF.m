% Run Stacking: broadening of response function due to different
% - retarding potentials
% - column densities
RunList = 'KNM2_Prompt';
TwinBias_Q = 18573.70;
savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_RunStackingRF_%s_E0%.2feV.mat',savedir,RunList,TwinBias_Q);

if exist(savename,'file')
    load(savename);
else
    RunAnaArg = {'RunList',RunList,... % all KNM2 golden runs
        'fixPar','E0 Bkg Norm',...           % free Parameter !!
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...
        'TwinBias_Q',TwinBias_Q,...
        'ROIFlag','14keV'};
    
    % build object of MultiRunAnalysis class
    A = MultiRunAnalysis(RunAnaArg{:});
    RFmodel = A.ModelObj.RF;
    
    RFsigma = std(A.SingleRunData.qU,0,2);
    A.ModelObj.MACE_Sigma = RFsigma;
    A.ModelObj.InitializeRF('RebinMode','Integral');
    RFmodelImp = A.ModelObj.RF;
    
    A.ModelObj.MACE_Sigma = 2.*RFsigma; % enhance
    A.ModelObj.InitializeRF('RebinMode','Integral');
    RFmodelImp2 = A.ModelObj.RF;
    
    A.ModelObj.MACE_Sigma = 5.*RFsigma; % enhance
    A.ModelObj.InitializeRF('RebinMode','Integral');
    RFmodelImp5 = A.ModelObj.RF;
    
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
    save(savename,'Te','qU','RF','RFmodel','RFmodelImp','RFmodelImp2','RFmodelImp5','RFmodelImp10',...
                  'RFsigma','qUModel','TeModel','nRuns','RunList','A');
    
    % interpola1 RF, so that they all have same Te-qU
    RFinter= zeros(nRuns,numel(TeModel),numel(qUModel));
    for i=1:nRuns
        progressbar(i/nRuns);
        for q=1:numel(qUModel)
            RFinter(i,:,q) = squeeze(interp1(Te(i,:)-qU(i,q),RF(i,:,q),TeModel-qUModel(q),'linear','extrap'));
        end
    end
    
    RFinterMean = squeeze(mean(RFinter));
    RFinterStd = squeeze(std(RFinter));
    save(savename,'RFinter','RFinterMean','RFinterStd','-append')
    
end
%% plot and compare
pSingleObj = cell(nRuns,1);
qUi = 27;
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
% for i=1:nRuns
%     pSingleObj{i} = plot(Te(i,:)-qU(i,qUi),RF(i,:,qUi),'LineWidth',1,'Color',rgb('PowderBlue'));
%     hold on;
% end
[l1,a1] = boundedline(TeModel-qUModel(qUi),RFinterMean(:,qUi),RFinterStd(:,qUi));
hold on;
l1.Color = rgb('DodgerBlue'); l1.LineWidth = 2;
a1.FaceColor = rgb('PowderBlue');
pref = plot(TeModel-qUModel(qUi),RFmodel(:,qUi),'-','Color',rgb('Orange'),'LineWidth',2);
pImp = plot(TeModel-qUModel(qUi),RFmodelImp(:,qUi),'--','Color',rgb('Crimson'),'LineWidth',2);
pImp2 = plot(TeModel-qUModel(qUi),RFmodelImp2(:,qUi),'-.','Color',rgb('DarkSlateGray'),'LineWidth',2);
pImp5 = plot(TeModel-qUModel(qUi),RFmodelImp5(:,qUi),'-.','Color',rgb('Magenta'),'LineWidth',2);
pImp10 = plot(TeModel-qUModel(qUi),RFmodelImp10(:,qUi),':','Color',rgb('YellowGreen'),'LineWidth',2);
xlabel(sprintf('{\\itE} - %.0f (eV)',qUModel(qUi)));
ylabel('Transmission probability');
PrettyFigureFormat('FontSize',24)
leg = legend([l1,a1,pref,pImp,pImp2,pImp5,pImp10],sprintf('Averaged RF     \\langleRF(x)\\rangle'),...
    sprintf('1\\sigma error band   \\langleRF(x)\\rangle'),...
    sprintf('Model average  RF(\\langlex\\rangle)'),...
    sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV',RFsigma(qUi)*1e3),...
    sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV',2*RFsigma(qUi)*1e3),...
    sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV',5*RFsigma(qUi)*1e3),...
    sprintf('RF(\\langlex\\rangle) with \\sigma = %.0f meV',10*RFsigma(qUi)*1e3));
leg.EdgeColor = rgb('Silver');
leg.Location = 'southeast';
xlim([2.78,3])





