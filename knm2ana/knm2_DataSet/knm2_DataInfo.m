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
Ba_Std = std(mean(d.A.SingleRunData.MACE_Ba_T,2));
qU = squeeze(d.A.SingleRunData.qU(1,1,:));
qU_std= std(qU);
