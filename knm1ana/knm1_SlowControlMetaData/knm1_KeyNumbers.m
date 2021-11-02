% some key parameters for knm1
% for PhD thesis
savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    return
end

if isempty(R.FitResult)
    R.Fit;
    save(savefile,'R','-append');
end
%% time
Sec2Hour = 60*60;
TimeAll = R.RunData.TimeSec/Sec2Hour;
Time90  = R.RunData.TimeSec*sum(R.RunData.qUfrac)/Sec2Hour;
Time40  = R.RunData.TimeSec*sum(R.RunData.qUfrac(R.exclDataStart:end))/Sec2Hour;
TimeBkg = R.RunData.TimeSec*sum(R.RunData.qUfrac(end-4:end))/Sec2Hour;

fprintf('Net time [-40,+50] eV : %.1f hours \n',Time40);
fprintf('Net time [-90,+50] eV : %.1f hours \n',Time90);
fprintf('Net time [-200,+50] eV: %.1f hours \n',TimeAll);
fprintf('Net time background %.1f hours (%.1f%% (90 eV), %.1f%% (40 eV))\n',TimeBkg,1e2*TimeBkg/Time90,1e2*TimeBkg/Time40);

%% number of electrons: analysis interval
BkgRate = R.ModelObj.BKG_RateSec;
Nbkg_qU  = R.RunData.qUfrac(R.exclDataStart:end).*R.RunData.TimeSec.*BkgRate;
Nbkg     = sum(Nbkg_qU);
Nall    = sum(R.RunData.TBDIS(R.exclDataStart:end));
Nsig = Nall-Nbkg;

fprintf('---------- 40 eV interval ----------------\n');
fprintf('Total      number electrons %.2e \n',Nall);
fprintf('Signal     number electrons %.2e \n',Nsig);
fprintf('Background number electrons %.2e \n',Nbkg);

%% signal to background
Nsig_qU = R.RunData.TBDIS(R.exclDataStart:end)-Nbkg_qU;
SB_qU = Nsig_qU./Nbkg_qU;
fprintf('--------------------------------\n');
fprintf('Tot signal to background   %.2f \n',Nsig/Nbkg);
fprintf('Mean signal to background   %.2f \n',mean(SB_qU));
fprintf('Min signal to background   %.2f \n',min(SB_qU));
fprintf('Max signal to background   %.2f \n',max(SB_qU));

