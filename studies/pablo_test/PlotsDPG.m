

addpath(genpath('../../../Samak2.0'));

opt_bkg = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',0.2};

A = TBD('nTeBinningFactor',100,'TimeSec',60*60,opt_bkg{:});
A.mnuSq_i = 0;



A.ComputeTBDDS();
A.ComputeTBDIS();

Am = TBD('nTeBinningFactor',100);
Am.mnuSq_i = (1).^2;

Am.ComputeTBDDS();

%Am.ComputeTBDIS();
figure(54669)
hold on
plot(A.qU-18575,A.TBDIS,'k');
bar(A.qU-18575,A.TBDIS);
   % Am.Te,Am.TBDDS);
   hold off
title('Tritium \beta-decay Spectrum 1 hour')
ylabel('count rate [a.u.]');
xlabel('E - E_0 [eV]');


A.AddStatFluctTBDIS();
figure(77)
hold on
bar(A.qU,A.TBDIS)
errorbar(A.qU,A.TBDIS,A.TBDISE*10,'k')
hold off
title('Tritium \beta-decay Spectrum 1 day with fluct.')
ylabel('count rate [a.u.]');
xlabel('E - E_0 [eV]');
legend({'','error bars \times 10'})


% 
% figure(2754)
% plot(A.Te(end-800:end-400)-18575,A.TBDDS(end-800:end-400),...
%     Am.Te(end-800:end-400)-18575,Am.TBDDS(end-800:end-400));
% ylabel('count rate [a.u.]');
% xlabel('E - E_0 [eV]');
