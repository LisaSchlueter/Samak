  %% configure RunAnalysis object
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType',DataType,...
        'fixPar','E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2_0p5eV',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope',3*1e-06,...
        'TwinBias_BKG_PtSlope',3*1e-06,...
        'DopplerEffectFlag','FSD'};
    A = MultiRunAnalysis(RunAnaArg{:});
    
    
    %%
    Time = A.ModelObj.qUfrac(11:end-5).*A.ModelObj.TimeSec;
    qU = abs(A.ModelObj.qU(11:end-5)-18574);
    GetFigure
    m4Sq = logspace(0,3,1e3);
    
   TBDIS = A.ModelObj.TBDIS(11:end-5);
   TBDIS_inter = interp1(qU.^2,TBDIS,linspace(min(m4Sq),max(m4Sq),1e4));
   
   TBDISE = sqrt(A.ModelObj.TBDIS(11:end-5));
  %  plot(qU.^2,1./TBDISE.^2,'x')
 %  plot(qU.^2,2./(2+TBDISE.^2),'x')

%set(gca,'YScale','log');
   % xlim([1 40])
    
   SigmaSqRel_Stat = 1./TBDISE;
   SigmaSqAbs_Stat = SigmaSqRel_Stat.*TBDIS;  

    SigmaSqRel_Syst_const = 0.0028;
    SigmaSqRel_Syst_lin = TBDIS.*0.05/max(TBDIS);%qU*0.0005;%1e-06;      

   SigmaSqAbs_Syst_const = SigmaSqRel_Syst_const.*TBDIS;
    SigmaSqAbs_Syst_lin = SigmaSqRel_Syst_lin.*TBDIS;
   
   SigmaSqAbs_Tot_const = SigmaSqAbs_Stat+SigmaSqAbs_Syst_const;
    SigmaSqAbs_Tot_lin = SigmaSqAbs_Stat+SigmaSqAbs_Syst_lin;
      
  p1 =  plot(qU.^2,SigmaSqAbs_Syst_const./SigmaSqAbs_Tot_const);
  hold on;
    p2 =  plot(qU.^2,SigmaSqAbs_Syst_lin./SigmaSqAbs_Tot_lin);

   legend('const','lin');
       set(gca,'XScale','log');
   
     xlim([3 1600])