function CheckmnuSqDist(varargin)
    p=inputParser;
    p.addParameter('mnuSq',0,@(x)isfloat(x));
    p.addParameter('pullFlag',3);
    p.addParameter('Nfit',1,@(x)isfloat(x));
    p.addParameter('DataType','Twin',@(x)ismember(x,{'Twin','Real'}));
    p.addParameter('Syst','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('RunList','KNM1',@(x)ischar(x));
    p.addParameter('Plot','ON',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    mnuSq    = p.Results.mnuSq;
    pullFlag = p.Results.pullFlag;
    Nfit     = p.Results.Nfit;
    DataType = p.Results.DataType;
    Syst     = p.Results.Syst;
    RunList  = p.Results.RunList;

    
    if exist(sprintf('./mNuFitResult_AllParams_mnuSq%g_Nfit%g.mat',mnuSq,Nfit),'file') && Nfit>1
        load(sprintf('./mNuFitResult_AllParams_mnuSq%g_Nfit%g.mat',mnuSq,Nfit));
    else
        if strcmp(Syst,'OFF')
            Chi2opt='chi2Stat';
            NP=1;
        else
            Chi2opt='chi2CMShape';
            NP=1.064;
        end
        
        fitresults = zeros(9,Nfit);
        
        for j=1:Nfit
            M = MultiRunAnalysis('RunList',RunList,...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
                        'chi2',Chi2opt,...              % uncertainties: statistical or stat + systematic uncertainties
                        'DataType',DataType,...                 % can be 'Real' or 'Twin' -> Monte Carlo
                        'fixPar','mNu E0 Norm Bkg',...    % free Parameter!!
                        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
                        'NonPoissonScaleFactor',NP,...     % background uncertainty are enhanced
                        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
                        'pullFlag',pullFlag,...
                        'FSDFlag','SibilleFull',...           % final state distribution
                        'ELossFlag','KatrinT2',...            % energy loss function
                        'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
                        'DopplerEffectFlag','OFF',...
                        'Twin_SameCDFlag','OFF',...
                        'Twin_SameIsotopFlag','OFF',...
                        'SynchrotronFlag','ON',...
                        'AngularTFFlag','OFF',...
                        'TwinBias_Q',18573.73,...
                        'TwinBias_mnuSq',mnuSq);

            if strcmp(DataType,'Twin')
                %M.RunData.TBDIS(end-5:end)=mean(M.RunData.TBDIS(end-5:end)./(M.RunData.qUfrac(end-5:end).*M.RunData.TimeSec)).*M.RunData.qUfrac(end-5:end).*M.RunData.TimeSec;
                %statfluct = zeros(numel(D.RunData.qU),1);
                %for i=1:numel(D.RunData.qU)
                %    gm=gmdistribution(D.RunData.TBDIS(i),(NP)*D.RunData.TBDIS(i));
                %    statfluct(i) = random(gm)-D.RunData.TBDIS(i);
                %end
                %M.RunData.TBDIS = M.RunData.TBDIS+statfluct;
                M.RunData.TBDIS = mvnrnd(M.RunData.TBDIS',M.FitCM,1)';
            end
            
            M.exclDataStart = M.GetexclDataStart(40);
            %M.InitModelObj_Norm_BKG('RecomputeFlag','ON');
            %M.ModelObj.BKG_RateSec_i=0.292256;
            M.Fit;
            fitresults(1,j)= M.FitResult.par(1);
            fitresults(2,j)= M.FitResult.err(1);
            fitresults(3,j)= M.ModelObj.Q_i+M.FitResult.par(2);
            fitresults(4,j)= M.FitResult.err(2);
            fitresults(5,j)= M.ModelObj.BKG_RateSec_i+M.FitResult.par(3);
            fitresults(6,j)= M.FitResult.err(3);
            fitresults(7,j)= M.FitResult.par(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1) + 1;
            fitresults(8,j)= M.FitResult.err(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1);
            fitresults(9,j)= M.FitResult.chi2min;
            save(sprintf('./mNuFitResult_AllParams_mnuSq%g_Nfit%g.mat',mnuSq,Nfit),'fitresults');
        end
    end