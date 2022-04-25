% do combi fit knm1+2
RecomputeFlag = 'OFF';
Plot = 'OFF';
tStart = tic;
%% KNM-1 Model
chi2          = 'chi2CMShape';
DataType      = 'Real';
KNM1SysBudget = 24;
KNM1Doppler   = 'OFF';
savedir = [getenv('SamakPath'),'knm12Combi/knm2_Combination/results/'];
savename = sprintf('%sknm2_CombiFit_%s_Uniform_%s_Knm1SysBudget%.0f_Knm1DE%s.mat',savedir,DataType,chi2,KNM1SysBudget,KNM1Doppler);

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    d = importdata(savename);
    fprintf('load results from file %s \n',savename)
else
    range                 = 40;
    freePar               = 'mNu E0 Bkg Norm';
    RunList = 'KNM1'; 
    K1 = MultiRunAnalysis('RunList',RunList,...
        'chi2',chi2,...
        'DataType',DataType,...
        'fixPar',freePar,... free parameter
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'SysBudget',KNM1SysBudget,...
        'ELossFlag','KatrinT2',...
        'SynchrotronFlag','ON',...
        'RadiativeFlag','ON',...
        'AngularTFFlag','OFF',...
        'DopplerEffectFlag',KNM1Doppler);
    
    K1.exclDataStart = K1.GetexclDataStart(range);
    K1.i_qUOffset    = zeros(1,K1.ModelObj.nPixels);
    K1.i_mTSq        = zeros(1,K1.ModelObj.nPixels);

    %% KNM-2 Model
    AnaFlag               = 'StackPixel';
    RingMerge             = 'Full';%'None'; %only relevand when AnaFlag = Ring
    DopplerEffectFlag     = 'FSD';
    BKG_PtSlope           = 3*1e-06;
    TwinBias_BKG_PtSlope  = 3*1e-06;
    FSDFlag               = 'KNM2_0p1eV';
    PullFlag              = 99;% 99 = no pull
    SysBudget             = 40;
    NonPoissonScaleFactor = 1.112;
    SigmaSq               =  0.0124+0.0025;
    
    % Init Model Object and covariance matrix object
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar',freePar,...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag',FSDFlag,...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',SysBudget,...
        'AnaFlag',AnaFlag,...
        'chi2',chi2,...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',NonPoissonScaleFactor,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge',RingMerge,...
        'PullFlag',PullFlag,...;%99 = no pull
        'BKG_PtSlope',BKG_PtSlope,...
        'TwinBias_BKG_PtSlope',TwinBias_BKG_PtSlope,...
        'DopplerEffectFlag',DopplerEffectFlag};
    K2 = MultiRunAnalysis(RunAnaArg{:});
    
    K2.exclDataStart = K2.GetexclDataStart(range);
    K2.i_qUOffset = zeros(1,K2.ModelObj.nPixels);
    K2.i_mTSq = zeros(1,K2.ModelObj.nPixels);
    
    
    %% Combine
    
    Data = [K1.RunData.qU(2:end),K2.RunData.qU,...
        K1.RunData.TBDIS(2:end),K2.RunData.TBDIS,...
        K1.RunData.TBDISE(2:end),K2.RunData.TBDISE];
    
    %%
  
    if strcmp(K1.chi2,'chi2Stat')
        if ~contains(K1.fixPar,'fix 3 ;')
            K1.InitModelObj_Norm_BKG('RecomputeFlag','OFF');
        end
        
        [StatCM, StatCMFrac] = K1.ComputeCM_StatPNP(varargin);
        K1.FitCM = StatCM;
        K1.FitCMShape = StatCM;
        K1.FitCMFrac = StatCMFrac;
        K1.FitCMFracShape = StatCMFrac;
        
        if ~contains(K2.fixPar,'fix 3 ;')
            K2.InitModelObj_Norm_BKG('RecomputeFlag','OFF');
        end
        
        [StatCM2, StatCMFrac2] = K2.ComputeCM_StatPNP(varargin);
        K2.FitCM      = StatCM2;
        K2.FitCMShape = StatCM2;
        K2.FitCMFrac = StatCMFrac2;
        K2.FitCMFracShape = StatCMFrac2;
    else
        K1.InitModelObj_Norm_BKG('RecomputeFlag','ON');
        K1.ComputeCM;
        K2.InitModelObj_Norm_BKG('RecomputeFlag','ON');
        K2.ComputeCM;
    end
    
    %% cov mats
    nqU = size(Data,1);
    COVMAT          = diag(zeros(2*nqU,2*nqU));
    COVMATShape     = diag(zeros(2*nqU,2*nqU));
    COVMATFracShape = diag(zeros(2*nqU,2*nqU));
    
    COVMAT(1:nqU,1:nqU) = K1.FitCM(2:end,2:end);
    COVMAT(nqU+1:2*nqU,nqU+1:2*nqU) = K2.FitCM;
    COVMATShape(1:nqU,1:nqU) = K1.FitCMShape(2:end,2:end);
    COVMATShape(nqU+1:2*nqU,nqU+1:2*nqU) = K2.FitCMShape;
    COVMATFracShape(1:nqU,1:nqU) = K1.FitCMFracShape(2:end,2:end);
    COVMATFracShape(nqU+1:2*nqU,nqU+1:2*nqU) = K2.FitCMFracShape;
    
    %%
    K1.ModelObj.nPixels = 2; % to make FITC read data properly
    
    F = FITC('SO',K1.ModelObj,...
        'SO2',K2.ModelObj,...
        'DATA',Data,'fitter',K1.fitter,...
        'chi2name',K1.chi2,'minuitOpt',K1.minuitOpt,...
        'COVMAT', real(COVMAT),'COVMATFrac', real(K1.FitCMFrac),...
        'COVMATShape', real(COVMATShape),'COVMATNorm',K1.FitCMNorm,...
        'COVMATFracShape',real(COVMATFracShape),...
        'pulls',K1.pulls,...
        'pullFlag',K1.pullFlag,...
        'fixPar','',...
        'exclDataStart',K1.exclDataStart,...
        'exclDataStop',K1.exclDataStop,...
        'i_mnu',K1.i_mnu,...
        'i_Q',K1.i_Q,...
        'i_Q2',K2.i_Q,...
        'i_B',K1.i_B,...
        'i_N',K1.i_N,...
        'i_DTGS',K1.i_DTGS,...
        'i_DTES',K1.i_DTES,...
        'i_HTGS',K1.i_HTGS,...
        'i_HTES',K1.i_HTES,...
        'i_TTGS',K1.i_TTGS,...
        'i_TTES',K1.i_TTES,...
        'i_qUOffset',K1.i_qUOffset,...
        'i_mTSq',K1.i_mTSq,...
        'i_FracTm',0);
    
    K1.FitResult = struct(...
        'par',F.RESULTS{1},....
        'err',F.RESULTS{2},....
        'chi2min',F.RESULTS{3},...
        'errmat',F.RESULTS{4},...
        'dof',F.RESULTS{5});
    
    tCPU = toc(tStart);
    MakeDir(savedir);
    
    FitResult = K1.FitResult;
    save(savename,'FitResult','K1','K2','tCPU','F')
    d = importdata(savename);
end


%% asym err
ErrNeg  = d.F.RESULTS{end-1};
ErrPos  = d.F.RESULTS{end};
MeanErr = 0.5.*(ErrPos-ErrNeg);
%% display
fprintf('- Fit results  --------------------\n');
fprintf('mnu^2      = %.2f +- %.2f eV^2 \n',d.FitResult.par(1),MeanErr(1));
fprintf('E0 (KNM-1) = %.2f +- %.2f eV \n',d.FitResult.par(2)+d.K1.ModelObj.Q_i,d.FitResult.err(2));
fprintf('E0 (KNM-2) = %.2f +- %.2f eV \n',d.FitResult.par(3)+d.K2.ModelObj.Q_i,d.FitResult.err(3));
fprintf('B  (KNM-1) = %.1f +- %.1f mcps \n',1e3.*(d.FitResult.par(4)+d.K1.ModelObj.BKG_RateSec_i),1e3*d.FitResult.err(4));
fprintf('B  (KNM-2) = %.1f +- %.1f mcps \n',1e3*(d.FitResult.par(5)+d.K2.ModelObj.BKG_RateSec_i),1e3*d.FitResult.err(5));
fprintf('N  (KNM-1) = %.3f +- %.3f  \n',(d.FitResult.par(6)+1),d.FitResult.err(6));
fprintf('N  (KNM-2) = %.3f +- %.3f  \n',(d.FitResult.par(7)+1),d.FitResult.err(7));
fprintf('chi^2_min  = %.1f (%.0f dof)\n',d.FitResult.chi2min,27+28-7)
fprintf('- ---------------------------------\n');

%% plot spectrum

if strcmp(Plot,'ON')
    FontSize = 22;
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.65]);
    s1 = subplot(4,1,1:3);
    Time1 = d.K1.ModelObj.TimeSec.*d.K1.ModelObj.qUfrac;
    data1 = errorbar(d.K1.RunData.qU-18574,d.K1.RunData.TBDIS./Time1,50.*sqrt(d.K1.RunData.TBDIS)./Time1,'.','LineWidth',2,'MarkerSize',15,'CapSize',0);
    hold on;
    fit1 = plot(d.K1.RunData.qU-18574,d.K1.ModelObj.TBDIS./Time1,'-','LineWidth',2,'Color',data1.Color);
    
    Time2 = d.K2.ModelObj.TimeSec.*d.K2.ModelObj.qUfrac;
    data2 = errorbar(d.K2.RunData.qU-18574,d.K2.RunData.TBDIS./Time2,50.*sqrt(d.K2.RunData.TBDIS)./Time2,'.','LineWidth',2,'MarkerSize',15,'CapSize',0);
    fit2 = plot(d.K2.RunData.qU-18574,d.K2.ModelObj.TBDIS./Time2,'-','LineWidth',2,'Color',data2.Color);
    
    PrettyFigureFormat('FontSize',FontSize);
    xlabel('Retarding potential - 18574 (eV)');
    ylabel('Rate (cps)');
    
    set(gca,'YScale','log');
    leg = legend([data1,data2,fit1,fit2],...
        sprintf('KNM1 data with 1\\sigma error bars \\times 50'),sprintf('KNM2 data with 1\\sigma error bars \\times 50'),...
        'Best-fit model','Best-fit model');
    PrettyLegendFormat(leg);
    leg.NumColumns = 2;
    ax1 = gca;
    Pos_i = ax1.Position;
    ax1.Position = [Pos_i(1),Pos_i(2)+0.05,Pos_i(3),Pos_i(4)];
    ylim([0.11,1e2]);
   xticklabels('');  
  
    s2 = subplot(4,1,4);
    pr1 = plot(d.K1.RunData.qU-18574,(d.K1.ModelObj.TBDIS-d.K1.RunData.TBDIS)./sqrt(d.K1.RunData.TBDIS),'.:','MarkerSize',15,'LineWidth',2,'Color',data1.Color);
    hold on;
    pr2 = plot(d.K2.RunData.qU-18574,(d.K2.ModelObj.TBDIS-d.K2.RunData.TBDIS)./sqrt(d.K2.RunData.TBDIS),'.:','MarkerSize',15,'LineWidth',2,'Color',data2.Color);
    leg = legend([pr1,pr2],'KNM-1','KNM-2');
    PrettyLegendFormat(leg);
    linkaxes([s1,s2],'x');
    xlim([-41,137])
    PrettyFigureFormat('FontSize',FontSize);
    xlabel('Retarding potential - 18574 (eV)');
    ylabel(sprintf('Residuals (\\sigma)'));
    ylim([-3 3]);
    ax2 = gca;
    Pos_i = ax2.Position;
    ax2.Position = [Pos_i(1),Pos_i(2)+0.01,Pos_i(3),Pos_i(4)+0.09];
    ax2.YLabel.Position(1) = ax1.YLabel.Position(1);
    %% save plot
    plotdir = strrep(savedir,'results','plots');
    plotname = sprintf('%sknm2_CombiFitUniform.pdf',plotdir);
    export_fig(gcf,plotname);%,'-dpng','-r300');
end