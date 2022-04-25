% calculate qU-Scan expectation with randomized MC
chi2 = 'chi2CMShape';
qURange  = [95,20];
NonPoissonScaleFactor = 1.112;
fitter = 'minuit';
nSamples = 1000;

savedir = [getenv('SamakPath'),'knm2ana/knm2_qUScan/results/'];
savename_randMC = sprintf('%sknm2_qUScan_MC_%.0feV_to_%.0feV_%s_NP%.3f_%s_%.0fsamples.mat',...
    savedir,qURange(1),qURange(2),chi2,NonPoissonScaleFactor,fitter,nSamples);

if exist(savename_randMC,'file')
    load(savename_randMC);
else
    SigmaSq =  0.0124+0.0025;
    BKG_PtSlope = 3*1e-06;
    ELossFlag = 'KatrinT2A20';
    FSDFlag = 'KNM2_0p1eV';
    SysBudget = 40;
    
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType','Twin',...
        'fixPar','mNu E0 Bkg Norm',...
        'RadiativeFlag','ON',...
        'DopplerEffectFlag','FSD',...
        'fitter',fitter,...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag',ELossFlag,...
        'SysBudget',SysBudget,...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_BKG_PtSlope',BKG_PtSlope};
    T = MultiRunAnalysis(RunAnaArg{:});
    range = 40;               % fit range in eV below endpoint
    T.exclDataStart = T.GetexclDataStart(range); % find correct data, where to cut spectrum.Fit;
    TBDIS_i = T.ModelObj.TBDIS;
    FitCM = T.FitCMShape;
    
    % randomize twin
    qU = T.ModelObj.qU;
    Time = T.ModelObj.qUfrac.*T.ModelObj.TimeSec;
    TBDIS_r = mvnrnd(TBDIS_i,FitCM,nSamples);
    exclDataStart_v = T.GetexclDataStart(qURange(1)):T.GetexclDataStart(qURange(2));
    nFits = numel(exclDataStart_v);
    
    %%
    D = copy(repmat(T,nSamples,1));
    D = reshape(D,numel(D),1);
    
    par = zeros(4,nFits,nSamples);
    err = zeros(4,nFits,nSamples);
    chi2min = zeros(nFits,nSamples);
    dof = zeros(nFits,1);
    
    ForFactor = 50;
    
    for k=1:ForFactor
        nSamples_Start = (k-1)*nSamples/ForFactor+1;
        nSamples_Stop = k*nSamples/ForFactor;
        nSamples_k = nSamples_Start:1:nSamples_Stop;
        
        savename_tmp = strrep(savename_randMC,'.mat',sprintf('_%.0fto%.0f.mat',nSamples_k(1),nSamples_k(end)));
        if exist(savename_tmp,'file')
            d = importdata(savename_tmp);
            par(:,:,nSamples_k)   = d.par(:,:,nSamples_k);
            err(:,:,nSamples_k)   = d.err(:,:,nSamples_k);
            chi2min(:,nSamples_k) = d.chi2min(:,nSamples_k);
        elseif ~exist(savename_tmp,'file')
            parfor i=nSamples_k
                D(i).SimulateStackRuns;
                D(i).RunData.TBDIS = TBDIS_r(i,:)';
                %D(i).Fit;
                
                [parqU, errqU, chi2qU, dofqU] =  D(i).qUScan(...
                    'qURange',qURange,...
                    'saveplot','OFF',...
                    'RecomputeFlag','ON',...
                    'CorrMean','OFF',...
                    'HoldOn','OFF',...
                    'RelFlag','OFF',...
                    'saveStr','',...
                    'RefLine','OFF',...
                    'saveResult','OFF',...
                    'PlotResult','OFF');
                
                par(:,:,i)   = parqU(1:4,:);
                err(:,:,i)   = errqU(1:4,:);
                chi2min(:,i) = chi2qU;
                dof = dofqU;
            end
            
            savename_tmp = strrep(savename_randMC,'.mat',sprintf('_%.0fto%.0f.mat',nSamples_k(1),nSamples_k(end)));
            save(savename_tmp,'par','err','chi2min','dof','nSamples_k');
        end
    end
    
    % save
    save(savename_randMC,'par','err','chi2min','dof',...
        'qU','Time','TBDIS_i','TBDIS_r',...
        'FitCM','exclDataStart_v','RunAnaArg','T','qURange');
    fprintf('save file to %s\n',savename_randMC);
end



