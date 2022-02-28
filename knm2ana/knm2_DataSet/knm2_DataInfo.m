% general information on knm2 data set

% get Data
savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2ubfinal_Fit_Bpng-%.1fmucpsPers_%s_%.0feV_%s_%s_%s_%s_SysBudget40.mat',...
    savedir,3,'Real',40,'mNuE0BkgNorm','chi2CMShape','StackPixel','KNM2_0p1eV');
d = importdata(savename);
%% number of electrons
CountsData_All = sum(d.A.RunData.TBDIS(d.A.exclDataStart:end));
CountsModel_All = sum(d.A.ModelObj.TBDIS(d.A.exclDataStart:end));
CountsModel_Bkg = sum(d.A.ModelObj.BKG_RateSec.*d.A.ModelObj.qUfrac(d.A.exclDataStart:end).*d.A.ModelObj.TimeSec);
CountsModel_Bkg_belowE0 = sum(d.A.ModelObj.BKG_RateSec.*d.A.ModelObj.qUfrac(d.A.exclDataStart:end-5).*d.A.ModelObj.TimeSec);
CountsModel_Signal = CountsModel_All-CountsModel_Bkg;

%% analyzing plane: magnetic field 
%d.A.ReadSingleRunData(); 
% Ba_Std = std(mean(d.A.SingleRunData.MACE_Ba_T,2));
% qU = squeeze(d.A.SingleRunData.qU(1,1,:));
% qU_std= std(qU);

%% display
% interval: 
qU = d.A.RunData.qU-18574;
qURM = d.A.RunData.qU_RM-18574;
fprintf('Measurement interval: qU-E0 = [%.0f ,%.0f] eV \n',qU(1),qU(end));
fprintf('Analysis interval:    qU-E0 = [%.0f ,%.0f] eV \n',qU(d.A.exclDataStart),qU(end));
fprintf('Rate monitor at:      qU-E0 = %.0f eV \n',qURM);

%% time
Sec2Hour = 60*60;
TimeAll = d.A.RunData.TimeSec/Sec2Hour;
Time90  = d.A.RunData.TimeSec*sum(d.A.RunData.qUfrac)/Sec2Hour;
Time40  = d.A.RunData.TimeSec*sum(d.A.RunData.qUfrac(d.A.exclDataStart:end))/Sec2Hour;
TimeBkg = d.A.RunData.TimeSec*sum(d.A.RunData.qUfrac(end-4:end))/Sec2Hour;

fprintf('------------------------ KNM2 ------------------------\n');
fprintf('%.0f golden scans \n',numel(d.A.SingleRunData.StartTimeStamp));
fprintf('First golden scan: %s \n',d.A.SingleRunData.StartTimeStamp(1));
fprintf('Last golden scan:  %s \n',d.A.SingleRunData.StartTimeStamp(end));
fprintf('Net time [-40,+135] eV : %.1f hours \n',Time40);
fprintf('Net time [-90,+135] eV : %.1f hours \n',Time90);
fprintf('Net time [-300,+135] eV: %.1f hours \n',TimeAll);
fprintf('Net time background %.1f hours (%.1f%% (90 eV), %.1f%% (40 eV))\n',TimeBkg,1e2*TimeBkg/Time90,1e2*TimeBkg/Time40);

%% number of electrons: analysis interval
BkgRate = d.A.ModelObj.BKG_RateSec;
Nbkg_qU  = d.A.RunData.qUfrac(d.A.exclDataStart:end).*d.A.RunData.TimeSec.*BkgRate;
Nbkg     = sum(Nbkg_qU);
Nall    = sum(d.A.RunData.TBDIS(d.A.exclDataStart:end));
Nsig = Nall-Nbkg;

fprintf('---------- 40 eV interval ----------------\n');
fprintf('Total      number electrons %.2e \n',Nall);
fprintf('Signal     number electrons %.2e \n',Nsig);
fprintf('Background number electrons %.2e \n',Nbkg);

%% signal to background
Nsig_qU = d.A.RunData.TBDIS(d.A.exclDataStart:end)-Nbkg_qU;
SB_qU = Nsig_qU./Nbkg_qU;
fprintf('--------------------------------\n');
fprintf('Tot signal to background   %.2f \n',Nsig/Nbkg);
fprintf('Mean signal to background   %.2f \n',mean(SB_qU));
fprintf('Min signal to background   %.2f \n',min(SB_qU));
fprintf('Max signal to background   %.2f \n',max(SB_qU));

%% activity
fprintf('Colun density %.2f x e17 \n',d.A.RunData.WGTS_CD_MolPerCm2*1e-17);
fprintf('Atomic purity %.2f %% \n',1e2*d.A.ModelObj.WGTS_epsT);


std(d.A.SingleRunData.WGTS_MolFrac_TT)
Activity = 2.*d.A.RunData.WGTS_CD_MolPerCm2.*d.A.ModelObj.WGTS_epsT.*pi*4.5^2.*d.A.ModelObj.TdecayC;
fprintf('Activity %.2e \n',Activity)
fprintf('TT %.2f %% \n',1e2*d.A.ModelObj.WGTS_MolFrac_TT);
fprintf('HT %.2f %% \n',1e2*d.A.ModelObj.WGTS_MolFrac_HT);
fprintf('DT %.2f %% \n',1e2*d.A.ModelObj.WGTS_MolFrac_DT);
fprintf('B source %.2f T \n',d.A.ModelObj.WGTS_B_T);
fprintf('B ana    %.2e T \n',d.A.ModelObj.MACE_Ba_T);
fprintf('B max    %.2f T \n',d.A.ModelObj.MACE_Bmax_T);