% calculate qU-Scan expectation with randomized MC

nSamples = 1000;
chi2 = 'chi2CMShape';
NP = 1.064;
qURange  = [95,20];
savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];
savename_randMC = sprintf('%sknm1_qUScan_MC_%.0feV_to_%.0feV_%s_NP%2g_%.0fsamples.mat',...
    savedir,qURange(1),qURange(2),chi2,NP,nSamples);

if exist(savename_randMC,'file')
    load(savename_randMC);
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
    dof = zeros(nFits,1);
    
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
        dof = dofqU;
    end
    
     
    % save
    save(savename_randMC,'par','err','chi2min','dof',...
        'qU','Time','TBDIS_i','TBDIS_r',...
        'FitCM','exclDataStart_v','RunAnaArg','T','qURange');
    fprintf('save file to %s\n',savename_randMC);
end

if ~exist('mNuSq_mu','var')
    
    nqU = arrayfun(@(x) numel(qU(x:end)),exclDataStart_v)';
    dof = nqU-4;
    
   
    %% next step: prepare results for plot
    mNuSq_abs = squeeze(par(1,:,:));               %(absolute fit values)
    mNuSq_40eV = mNuSq_abs(exclDataStart_v(13),:); % result in 40 eV range (standard analysis interval);
    mNuSq_rel = mNuSq_abs-mNuSq_40eV;              % relative to 40eV result
    
    % remove outliers further than 10sigma (1 sigma taken from of smallest qU-range fit ~2.2 eV);
   % Idx_rm = (abs(mNuSq_rel(1:end,:))>22);
    Idx_rm = abs(mNuSq_abs)>22; 
    mNuSq_rel(Idx_rm) = NaN;
    fprintf('remove outliers: %.1f%% \n',1e2*sum(sum(Idx_rm))./numel(Idx_rm));
    
    mNuSq_mu = nanmean(mNuSq_rel,2);
    mNuSq_std =nanstd(mNuSq_rel,0,2);
    
    GetFigure;
    boundedline(qU(exclDataStart_v)-18574,mNuSq_mu,mNuSq_std);
    grid on
    
    E0_abs = squeeze(par(2,:,:))+T.ModelObj.Q_i;               %(absolute fit values)
    E0_40eV = E0_abs(exclDataStart_v(13),:); % result in 40 eV range (standard analysis interval);
    E0_rel = E0_abs-E0_40eV;              % relative to 40eV result
    % remove the same outliers
    E0_rel(Idx_rm) = NaN; 
    E0_mu = nanmean(E0_rel,2);
    E0_std =nanstd(E0_rel,0,2);
    
    B_abs = squeeze(par(3,:,:))+T.ModelObj.BKG_RateSec_i ;     %(absolute fit values)
    B_40eV = B_abs(exclDataStart_v(13),:); % result in 40 eV range (standard analysis interval);
    B_rel = B_abs-B_40eV;              % relative to 40eV result
    % remove the same outliers
    B_rel(Idx_rm) = NaN; 
    B_mu = nanmean(B_rel,2);
    B_std =nanstd(B_rel,0,2);
    
    N_abs = squeeze(par(4,:,:))+1;     %(absolute fit values)
    N_40eV = N_abs(exclDataStart_v(13),:); % result in 40 eV range (standard analysis interval);
    N_rel = N_abs-N_40eV;              % relative to 40eV result
    % remove the same outliers
    N_rel(Idx_rm) = NaN; 
    N_mu = nanmean(N_rel,2);
    N_std =nanstd(N_rel,0,2);
    
    % p-value doesn't make that much sense, because not gaussian. will not include in display
    p_abs = 1-chi2cdf(chi2min,repmat(dof,1,nSamples));  
    p_abs(Idx_rm) = NaN;  % remove the same outliers
    p_mu = nanmean(p_abs,2);
    p_std =nanstd(p_abs,0,2);
    
    save(savename_randMC,...
        'mNuSq_mu','mNuSq_std','mNuSq_rel','mNuSq_abs',...
        'E0_mu','E0_std','E0_rel','E0_abs',...
        'B_mu','B_std','B_rel','B_abs',...
        'N_mu','N_std','N_rel','N_abs',...
        'dof','nqU',...
        'p_abs','p_mu','p_std',...
        'Idx_rm','-append');
end






