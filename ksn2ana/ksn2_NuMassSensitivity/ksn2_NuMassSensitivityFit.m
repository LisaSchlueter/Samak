% get nu-mass sensitivity with fits -> chi^2 profile


chi2 = 'chi2CMShape';
DataType = 'Real';%'Twin';
SavePlt = 'ON';

if strcmp(DataType,'Real')
    mNuSq = -1:0.2:2;%.5;
else
    mNuSq = -1:0.1:1;
end

PullFlag = 27;
mNu4Sq_all = [100,50,10,1];%logspace(-1,log10(40^2),nGridSteps);
sin2T4_all = [0.02,0.01,0.1,0.5];%logspace(-3,log10(0.5),nGridSteps);

for k = 1:numel(mNu4Sq_all)
    for j =1:numel(sin2T4_all)
        mNu4Sq_i = mNu4Sq_all(k); % 95;%logspace(-1,log10(40^2),nGridSteps);
        sin2T4_i = sin2T4_all(j);%0.02;%logspace(-3,log10(0.5),nGridSteps);
        
        
        savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuMassSensitivity/results/'];
        savefile = sprintf('%sksn2_NuMassSensitivityFits_%s_Pull%.0f_%s_mNuSqMin%.2f_mNuSqMax%.2f_mNuSteps%.0f.mat',...
            savedir,DataType,PullFlag,chi2,min(mNuSq),max(mNuSq),numel(mNuSq));
        
        if sin2T4_i~=0 || mNu4Sq_i~=0
            savefile = strrep(savefile,'.mat',sprintf('_sin2T4i%.3g_mNu4Sqi%.1feV2.mat',sin2T4_i,mNu4Sq_i));
        end
        
        if exist(savefile,'file')
            load(savefile,'mNuSq','chi2min');
        else
            
            %% configure RunAnalysis object
            if strcmp(chi2,'chi2Stat')
                NonPoissonScaleFactor = 1;
            elseif  strcmp(chi2,'chi2CMShape')
                NonPoissonScaleFactor = 1.112;
            end
            
            RunAnaArg = {'RunList','KNM2_Prompt',...
                'DataType',DataType,...
                'fixPar','E0 Norm Bkg sin2T4 mnu4Sq',...%free par
                'SysBudget',40,...
                'fitter','minuit',...
                'minuitOpt','min;migrad',...
                'RadiativeFlag','ON',...
                'FSDFlag','KNM2_0p5eV',...
                'ELossFlag','KatrinT2A20',...
                'AnaFlag','StackPixel',...
                'chi2',chi2,...
                'NonPoissonScaleFactor',NonPoissonScaleFactor,...
                'FSD_Sigma',sqrt(0.0124+0.0025),...
                'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                'TwinBias_Q',18573.7,...
                'PullFlag',PullFlag,...;%99 = no pull
                'BKG_PtSlope',3*1e-06,...
                'TwinBias_BKG_PtSlope',3*1e-06,...
                'DopplerEffectFlag','FSD'};
            A = MultiRunAnalysis(RunAnaArg{:});
            A.exclDataStart = A.GetexclDataStart(40);
            
            FitResults = cell(numel(mNuSq),1);
            chi2min    = zeros(numel(mNuSq),1);
            
            progressbar('Nu-mass sensitivity...');
            for i = 1:numel(mNuSq)
                progressbar((i-1)./numel(mNuSq));
                A.ModelObj.SetFitBias(0); % reset
                A.ModelObj.SetFitBiasSterile(mNu4Sq_i,sin2T4_i);
                A.ModelObj.mnuSq_i = mNuSq(i);
                A.Fit;
                FitResults{i} = A.FitResult;
                chi2min(i) = A.FitResult.chi2min;
            end
            
            save(savefile,'FitResults','mNuSq','chi2min','RunAnaArg','mNu4Sq_i','sin2T4_i');
            
        end
    end
end

%%
GetFigure;
plot(mNuSq,chi2min-min(chi2min),'LineWidth',2);

sin2T4 = cell2mat(cellfun(@(x) x.par(16),FitResults,'UniformOutput',0));
mNu4Sq =  cell2mat(cellfun(@(x) x.par(15),FitResults,'UniformOutput',0));

hold on;
plot(mNuSq,mNu4Sq-min(mNu4Sq),'LineWidth',2);
plot(mNuSq,sin2T4-min(sin2T4),'LineWidth',2);
xlabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
ylabel(sprintf('{\\itx} - {\\itx}_{min}'));
legend(sprintf('\\chi^2_{min}'),sprintf('{\\itm}_4^2 (eV^2)'),sprintf('|{\\itU}_{e4}|^2'),'Location','northwest');
legend boxoff
PrettyFigureFormat;