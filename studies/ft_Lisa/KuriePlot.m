function KuriePlot(varargin)
addpath(genpath('../../../Samak2.0'));

p = inputParser;
p.addParameter('TBDObj','', @(x) isa(x,'TBD'));
p.addParameter('Bias',0,@(x)isfloat(x));
p.addParameter('Data',0,@(x)isfloat(x));
p.parse(varargin{:}); 

TBDObj = p.Results.TBDObj;
Bias   = p.Results.Bias;
Data   = p.Results.Data;
% Model Fit
TBDObj.ComputeTBDDS(...
                'mSq_bias',Bias(1),...
                'E0_bias',Bias(2),...
                'B_bias',Bias(3),...
                'N_bias',Bias(4));
TBDObj.ComputeTBDIS;
TBDISFit = TBDObj.TBDIS;

% Model mNu= 0 eV
TBDObj.ComputeTBDDS(...
                'mSq_bias',0,...
                'E0_bias',0,...
                'B_bias',Bias(3),...
                'N_bias',Bias(4));
TBDObj.ComputeTBDIS;
TBDIS_ref = TBDObj.TBDIS;

% Soectrum - BGK (Rate)
BFitRate = (TBDObj.BKG_RateSec_i+Bias(3));
DataRate = Data(:,2)./(TBDObj.TimeSec.*TBDObj.qUfrac)-BFitRate; DataRate(DataRate<0)=0;
FitRate = TBDISFit./(TBDObj.TimeSec.*TBDObj.qUfrac)-BFitRate; FitRate(FitRate<0)=0;
RefRate =  TBDIS_ref./(TBDObj.TimeSec.*TBDObj.qUfrac)-BFitRate; RefRate(RefRate<0)=0;

FermiFunc = TBDObj.ComputeFermiCorr;
FermiFunc_qU = interp1(TBDObj.Te,FermiFunc,TBDObj.qU);
DataKurie = sqrt(DataRate./FermiFunc_qU);
FitKurie =  sqrt(FitRate./FermiFunc_qU);
RefKurie =  sqrt(RefRate./FermiFunc_qU);
%DataKurieErr = 0.5./sqrt(DataRate.*FermiFunc_qU).*sqrt(Data);

figure(999)
plot(Data(:,1),DataKurie,'x');
hold on;
plot(TBDObj.qU,FitKurie,'-')
plot(TBDObj.qU,RefKurie,'--')
hold off
%plot(TBDObj.qU,FitKurie,'-')

%plot(TBDObj.Te,FitKurie,'-')
%ylim([min(FitKurie) max(FitKurie)]);
%xlim([max(TBDObj.qU)-100 max(TBDObj.qU)])
end