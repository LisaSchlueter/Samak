% calculate qU-Scan expectation with randomized MC

nSamples = 3;%1000;
chi2 = 'chi2CMShape';
NP = 1.064;
qURange  = [95,90];%[95,20];
savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];
savename = sprintf('%sknm1_qUScan_MC_%.0feV_to_%.0feV_%s_NP%2g_%.0fsamples.mat',...
    savedir,qURange(1),qURange(2),chi2,NP,nSamples);

if exist(savename,'file')
else
    MakeDir(savedir);
    RunAnaArg = {'RunList','KNM1',...
        'RingMerge','None',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'DataType','Twin',...
        'fixPar','mNu E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',NP,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'SysBudget',22,...
        'AngularTFFlag','OFF'};
    
    T = MultiRunAnalysis(RunAnaArg{:});
    %
    T.exclDataStart = T.GetexclDataStart(40);
    T.Fit;
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
     
    parfor i=1:nSamples
        D(i).SimulateStackRuns;
        D(i).exclDataStart = D(i).GetexclDataStart(40);
        D(i).RunData.TBDIS = TBDIS_r(i,:)';
        D(i).FitCM = FitCM;
        D(i).Fit;
        
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
    end
    
    dof = dofqU;
     
    % save
    save(savename,'par','err','chi2min','dof',...
        'qU','Time','TBDIS_i','TBDIS_r',...
        'FitCM','exclDataStart_v','RunAnaArg','T','qURange');
    fprintf('save file to %s\n',savename);
end
