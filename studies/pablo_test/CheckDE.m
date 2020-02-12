clear; 
addpath(genpath('../../../Samak2.0'));


def = 'DopplerEffectFlag';
nte = 'nTeBinningFactor';
kk = 1; maxnte = 100;
for ii = 1:10:100
nten = ii;
A_DE = ref_ftmultipix(def,'matConv','TD','FT-TL3',...
    'KTFFlag','Compute','UseParallelRF','ON',...
    nte,nten);

A_OFF = ref_ftmultipix(def,'OFF','TD','Flat30',...
    'KTFFlag','Compute','UseParallelRF','ON',...
    nte,nten);

A_DE.ComputeTBDDS();
A_DE.ComputeTBDIS();

A_OFF.ComputeTBDDS();
A_OFF.ComputeTBDIS();

normcheck(kk) = A_DE.NormFactorTBDDS./A_OFF.NormFactorTBDDS;
TBDISratio = A_DE.TBDIS./(A_OFF.TBDIS*normcheck(kk));
ratiocheck(kk) = max(TBDISratio);
kk = kk + 1;
disp(kk)
end

plot(1:10:maxnte,ratiocheck)
%plot(A_DE.qU,TBDISratio)
