function ProfileLikelihood(varargin)
    p=inputParser;
    p.addParameter('RunList','KNM1',@(x)ismember(x,{'KNM1','KNM2_Prompt'}));
    p.addParameter('Parameter','eta',@(x)ismember(x,{'eta','mNu'}));
    p.addParameter('SysBudget',24,@(x)isfloat(x));
    p.addParameter('Nfit',100,@(x)isfloat(x));
    p.parse(varargin{:});
    RunList   = p.Results.RunList;
    Parameter = p.Results.Parameter;
    SysBudget = p.Results.SysBudget;
    Nfit      = p.Results.Nfit;

    matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Misc/')];
    savename=[matFilePath,sprintf('LikelihoodFct%s_%s_SysBudget%i_Nfit%i.mat',RunList,Parameter,SysBudget,Nfit)];
    if exist(savename,'file')
        load(savename);
        figure(1);
        plot(ScanRange,P,'LineWidth',2);
        cdf=zeros(numel(P),2);
        cdf(:,1) = ScanRange;
        for i=1:numel(P)
            cdf(i,2)=sum(P(1:i))./sum(P);
        end
        load('EtaFitResult_KNM2_Prompt_AllParams_mnuSq0_Nfit1000.mat','fitresults');
        fitresults(9,:)=fitresults(9,:)-1.6e10;
        [~,p]=kstest(fitresults(9,:),'CDF',cdf);
        figure(2);
        cdfplot(fitresults(9,:));
        hold on;
        x_values = linspace(min(fitresults(9,:)),max(fitresults(9,:)));
        plot(cdf(:,1),cdf(:,2),'r-');
        legend('Empirical CDF','Profile Likelihood','Location','best');
        title(sprintf('p-value: %.1g',p));
    else
        switch Parameter
            case 'mNu'
                fitPar='E0 Norm Bkg';
                ScanRange=linspace(-4,4,Nfit);
            case 'eta'
                fitPar='mNu E0 Norm Bkg';
                switch RunList
                    case 'KNM1'
                        ScanRange=linspace(-8e11,8e11,Nfit);
                    case 'KNM2_Prompt'
                        ScanRange=linspace(-4e11,4e11,Nfit);
                end
        end
        Chi2=ScanRange;

        if strcmp(RunList,'KNM1')
            D = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                'fixPar',fitPar,...                   % free Parameter!!
                'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
                'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                'pullFlag',99,...
                'FSDFlag','SibilleFull',...          % final state distribution
                'ELossFlag','KatrinT2',...            % energy loss function
                'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
                'DopplerEffectFlag','FSD',...
                'Twin_SameCDFlag','OFF',...
                'Twin_SameIsotopFlag','OFF',...
                'SynchrotronFlag','ON',...
                'AngularTFFlag','OFF',...
                'TwinBias_Q',18573.73,...
                'TwinBias_mnuSq',0);
        elseif strcmp(RunList,'KNM2_Prompt')
            D = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                'chi2','chi2CMShape',...                 % uncertainties: statistical or stat + systematic uncertainties
                'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
                'fixPar',fitPar,...                   % free Parameter!!
                'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                'NonPoissonScaleFactor',1.112,...     % background uncertainty are enhanced
                'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                'pullFlag',99,...
                'FSDFlag','KNM2',...          % final state distribution
                'ELossFlag','KatrinT2A20',...            % energy loss function
                'SysBudget',40,...                    % defines syst. uncertainties -> in GetSysErr.m;
                'DopplerEffectFlag','FSD',...
                'Twin_SameCDFlag','OFF',...
                'Twin_SameIsotopFlag','OFF',...
                'SynchrotronFlag','ON',...
                'AngularTFFlag','ON',...
                'TwinBias_Q',18573.7,...
                'TwinBias_mnuSq',0,...
                'FSD_Sigma',sqrt(0.0124+0.0025),...
                'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
                'BKG_PtSlope',3*1e-06,...
                'TwinBias_BKG_PtSlope',3*1e-06);
        end
        
        D.exclDataStart = D.GetexclDataStart(40);

        for i=1:numel(ScanRange)
            switch Parameter
                case 'mNu'
                    D.ModelObj.mnuSq_i=ScanRange(i);
                case 'eta'
                    D.ModelObj.eta_i = ScanRange(i);
            end
            D.Fit;
            Chi2(i)=D.FitResult.chi2min;
        end

        P = exp(-0.5.*Chi2);
        P = P./simpsons(ScanRange,P);
        save(savename,'ScanRange','P');
        cdf=zeros(numel(P),2);
        cdf(:,1) = ScanRange;
        for i=1:numel(P)
            cdf(i,2)=sum(P(1:i))./sum(P);
        end
        load('EtaFitResult_KNM2_Prompt_AllParams_mnuSq0_Nfit1000.mat','fitresults');
        %fitresults(9,:)=fitresults(9,:)-2.4e10;
        [~,p]=kstest(fitresults(9,:),'CDF',cdf);
        figure(2);
        cdfplot(fitresults(9,:));
        hold on;
        x_values = linspace(min(fitresults(9,:)),max(fitresults(9,:)));
        plot(cdf(:,1),cdf(:,2),'r-');
        legend('Empirical CDF','Profile Likelihood','Location','best');
        title(sprintf('p-value: %.1g',p));
    end