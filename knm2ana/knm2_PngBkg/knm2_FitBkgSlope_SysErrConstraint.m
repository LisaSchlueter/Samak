


BkgSlopeSigma_all = sort([0.5,1:2:20].*1e-06);%;4.74.*1e-06; %
DataType = 'Real';

savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2_FitBkgSlope_%s_SysErr_BkgConstraints_min%2g_max%2g_nfit%.0f.mat',...
    savedir,DataType,min(BkgSlopeSigma_all).*1e6,max(BkgSlopeSigma_all).*1e6,numel(BkgSlopeSigma_all));

if exist(savename,'file')
    d = importdata(savename);
else
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'chi2','chi2Stat',...
        'DataType',DataType,...
        'fixPar','mNu E0 Bkg Norm',...
        'RadiativeFlag','ON',...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',40,...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'RingMerge','Full',...
        'PullFlag',99,...
        'BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(40);
    
    if strcmp(DataType,'Twin')
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    %%
    savenameStat = sprintf('%sknm2_FitBkgSlope_%s_StatErr.mat', savedir,DataType);
    
    if exist(savenameStat,'file')
        load(savenameStat,'FitResultStat');
    else
        MakeDir(savedir);
        A.Fit;
        FitResultStat = A.FitResult;
        save(savenameStat,'FitResultStat');
    end
    
    mNuSqerrStat = 0.5*(FitResultStat.errPos(1)-FitResultStat.errNeg(1));
    mNuSqErr_all = zeros(numel(BkgSlopeSigma_all),1);
    mNuSqErr_all_sym = zeros(numel(BkgSlopeSigma_all),1);
    
    mNuSqStat = FitResultStat.par(1);
    mNuSq_all = zeros(numel(BkgSlopeSigma_all),1);
    FitResults = cell(numel(BkgSlopeSigma_all),1);
    %% fit with free bkg slope
    
    A.fixPar = ConvertFixPar('freePar','mNu E0 Bkg Norm BkgSlope','nPar',17);
    A.pullFlag = 10;% add background slope constraint with variable sigma
    progressbar('fitting bkg slope...');
    for i=1:numel(BkgSlopeSigma_all)
        progressbar(i/numel(BkgSlopeSigma_all));
        savename_tmp = sprintf('%sknm2_FitBkgSlope_%s_SysErr_BkgConstraints_%2g.mat',...
            savedir,DataType,BkgSlopeSigma_all(i).*1e6);
        if exist(savename_tmp,'file') 
            dtmp = importdata(savename_tmp);
            mNuSqErr_all(i) = dtmp.mNuSqErr;
            mNuSqErr_all_sym(i) = dtmp.FitResult.err(1);
            mNuSq_all(i) = dtmp.mNuSq;
            FitResults{i} = dtmp.FitResult;
        else
            A.ModelObj.SetFitBias(0)
            A.pulls =    BkgSlopeSigma_all(i);
            A.Fit;
            FitResult = A.FitResult;
            mNuSqErr = 0.5*(FitResult.errPos(1)-FitResult.errNeg(1));
            mNuSq = FitResult.par(1);
            save(savename_tmp,'FitResult','mNuSqErr','mNuSq')
            
            mNuSq_all(i) = mNuSq;
            mNuSqErr_all(i) = mNuSqErr;
            FitResults{i} = FitResult;
            mNuSqErr_all_sym(i) = FitResult.err(1);
        end
        
    end
    
    save(savename,'mNuSqerrStat','mNuSqStat','mNuSqErr_all','mNuSqErr_all_sym','mNuSq_all','BkgSlopeSigma_all','FitResults')
end

%% plot
%mNuSqErr_all = 
mNuSqerr_sys = sqrt(d.mNuSqErr_all.^2-d.mNuSqerrStat.^2);
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.4]);
x = linspace(0,max(d.BkgSlopeSigma_all),1e2);
y = interp1([0,d.BkgSlopeSigma_all],[0,mNuSqerr_sys'],x,'spline');

pref = plot(4.74*3.*ones(1e2,1),linspace(0,0.2,1e2),':','LineWidth',2,'Color',rgb('Gray'));
hold on;
plot(4.74.*ones(1e2,1),linspace(0,0.2,1e2),':','LineWidth',2,'Color',pref.Color);
pK2 = plot(linspace(0,30,1e2),0.16.*ones(1e2,1),'--','LineWidth',2,'Color',rgb('Orange'));
plot(x.*1e6,y,'-','LineWidth',3,'Color',rgb('DodgerBlue'));
xlabel(sprintf('\\sigma({\\its}_{qU}) (mcps/keV)'));
ylabel(sprintf('\\sigma_{syst.}({\\itm}_\\nu^2) (eV^{ 2})'));
PrettyFigureFormat;
xlim([0 15])
ylim([0 0.17]);
t1 = text(13.95,0.017,sprintf('3 \\times First Tritium'),'Rotation',90,'FontSize',get(gca,'FontSize'),'Color',pref.Color);
t1 = text(4.48,0.07,sprintf('First Tritium'),'Rotation',90,'FontSize',get(gca,'FontSize'),'Color',pref.Color);

tk2 = text(8,0.15,sprintf('{\\it s}_{qU} unconstrained'),'FontSize',get(gca,'FontSize'),'Color',pK2.Color);

pltdir = strrep(savedir,'results','plots');
MakeDir(pltdir);
pltname = [pltdir,'BkgSlopeSysErr.pdf'];
export_fig(pltname);



