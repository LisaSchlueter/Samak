addpath(genpath('../../../Samak2.0'));

myOpt= {'CPS','ON', 'FPD_Segmentation','PIXEL','FPD_Pixel',1,'TD','KrL3_32_HS','KTFFlag', 'MACER'...
       'DopplerEffectFlag','ON','ConvFlag','DopplerEffect', 'HVRipples','ON',...
       'L3_32_Phi0_i',100, 'BKG_RatePixelSec', 15}; %Activity and Offset to be modified
%Energy and Width to be modified in Kr.m class directly

HVRipplesInit = 0.52; %Initial value for HVRipple
global KrObjSim;
KrObjSim = InitKrKATRIN_Fit(myOpt{:},'HVRipplesP2PV', HVRipplesInit);
KrObjSim.ComputeKrDS; 
HVRipplesValue = [0.0001:0.0001:1];%Scan HVRipple Values in 0.1 mV steps

myKrIS      = zeros(1e4,KrObjSim.nqU);  %for storing IS and ISE
myKrISE     = zeros(1e4,KrObjSim.nqU);
tic

progressbar('Generate KrIS');
for i= 1:1e4 % (-0.5199*nSpec):(0.4800*nSpec)%-5199:4800 %10,000 spectra        
progressbar(i/1e4);
KrObjSim.HVRipplesP2PV = HVRipplesValue(i) %HVRipple: 0.1 - 1000 mV
KrObjSim.ComputeKrIS;
myKrIS(i,:) = KrObjSim.KrIS;
myKrISE(i,:) = KrObjSim.KrISE;
end
toc


MyCovMat = cov(myKrIS);
%save('../../studies/krypton_fit/CovMat_HVRippleSystematics.txt','MyCovMat','-ascii');

