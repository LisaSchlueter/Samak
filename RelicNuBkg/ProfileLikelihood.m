function ProfileLikelihood(varargin)
    p=inputParser;
    p.addParameter('Parameter','eta',@(x)ismember(x,{'eta','mNu'}));
    p.addParameter('SysBudget',24,@(x)isfloat(x));
    p.addParameter('Nfit',20,@(x)isfloat(x));
    p.parse(varargin{:});
    Parameter = p.Results.Parameter;
    SysBudget = p.Results.SysBudget;
    Nfit      = p.Results.Nfit;

    switch Parameter
        case 'mNu'
            fitPar='E0 Norm Bkg';
            ScanRange=linspace(-4,4,Nfit);
        case 'eta'
            fitPar='mNu E0 Norm Bkg';
            ScanRange=linspace(-8e11,8e11,Nfit);
    end
    Chi2=ScanRange;

    D = MultiRunAnalysis('RunList','KNM1',...       % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chi2CMShape',...                    % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Twin',...                       % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar',fitPar,...                         % free Parameter!!
        'RadiativeFlag','ON',...                    % theoretical radiative corrections applied in model
        'NonPoissonScaleFactor',1.064,...           % background uncertainty are enhanced
        'minuitOpt','min ; minos',...               % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...                 % final state distribution
        'ELossFlag','KatrinT2',...                  % energy loss function
        'SysBudget',SysBudget,...                   % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','OFF',...
        'Twin_SameCDFlag','OFF',...
        'Twin_SameIsotopFlag','OFF',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'TwinBias_Q',18573.73,...
        'TwinBias_mnuSq',0);
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
    matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Misc/')];
    savename=[matFilePath,sprintf('LikelihoodFctKRN1_%s_SysBudget%i_Nfit%i.mat',Parameter,SysBudget,Nfit)];
    save(savename,'ScanRange','P');