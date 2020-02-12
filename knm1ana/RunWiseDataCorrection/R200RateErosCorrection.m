function RunWiseErosFactor =  R200RateErosCorrection(varargin)
% Compute a Run Wise Correction Factor
% R200 is the reference rate at -200 eV below E0
% Cancel Correlations between R200 and
%  - Column Density
%  - [TT]
%  - [HT]
%  - [DT]
%  - qU200
%
% The corrected Rates is a linear combination written as:
% Rcorrected = R -
%              alpha . (rhoD - MeanrhoD)
%              beta  . (TT - MeanTT)
%              gamma . (HT - MeanHT)
%              eta   . (DT - MeanDT)
%              theta . (qU - MeanqU)
% 
% All subruns from a single run undergo the same correction
% Return a Run-wise Mutliplicative Correction factor
% 
% Method Adapted from EROS2 experiment, T. Lasserre PhD, 2000
% T. Lasserre, Last Modified 02/05/2019
% 

%% Bricolage pour récupérer R200/qU200
option = {...
    'DataType','Real',...
    'RunList','KNM1',...
    'exclDataStart',2,...
    'fixPar','1 5 6 7 8 9 10',...
    'chi2','chi2Stat',...
    'ELossFlag','KatrinD2',...
    'StackTolerance',1,...
    'NonPoissonScaleFactor',1.1,...
    'Debug','ON',...
    'AnaFlag','StackPixel',...
    'PixList',[]...
    };
SPR=MultiRunAnalysis(option{:});
R200    = (SPR.SingleRunData.TBDIS(1,:)./SPR.SingleRunData.TimeperSubRunperPixel(1,:,1))';
qU200   = SPR.SingleRunData.qU(1,:);

%% Definition of Parameters in use
p1=SPR.SingleRunData.WGTS_CD_MolPerCm2';
p2=SPR.SingleRunData.WGTS_MolFrac_TT';
p3=SPR.SingleRunData.WGTS_MolFrac_HT';
p4=SPR.SingleRunData.WGTS_MolFrac_DT';
p5=qU200';

m1=mean(p1);std1=std(p1);
m2=mean(p2);std2=std(p2);
m3=mean(p3);std3=std(p3);
m4=mean(p4);std4=std(p4);
m5=mean(p5);std5=std(p5);

rhoR200p1=corr(R200,p1); 
rhoR200p2=corr(R200,p2);  
rhoR200p3=corr(R200,p3);  
rhoR200p4=corr(R200,p4);  
rhoR200p5=corr(R200,p5);

rho12=corr(p1,p2);
rho13=corr(p1,p3);
rho23=corr(p2,p3);
rho14=corr(p1,p4);
rho24=corr(p2,p4);
rho34=corr(p3,p4);
rho15=corr(p1,p5);
rho25=corr(p2,p5);
rho35=corr(p3,p5);
rho45=corr(p4,p5);

%% Solving Equation for Cancelling Correlations Simultaneously
M = [...
      std1         rho12*std2     rho13*std3    rho14*std4  rho15*std5      ; ...
      rho12*std1   std2           rho23*std3    rho24*std4  rho25*std5      ; ...
      rho13*std1   rho23*std2     std3          rho34*std4  rho35*std5      ; ...
      rho14*std1   rho24*std2     rho34*std3    std4        rho45*std5      ; ...
      rho15*std1   rho25*std2     rho35*std3    rho45*std4  std5            ; ...
    ];
X = std(R200) * [rhoR200p1 rhoR200p2 rhoR200p3 rhoR200p4 rhoR200p5]';
S = M\X;

RunWiseErosFactor = (R200  - S(1) * (p1-m1) - S(2) * (p2-m2) - S(3) * (p3-m3) - S(4) * (p4-m4) - S(5) * (p5-m5))./R200;

end
