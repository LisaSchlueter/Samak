clear;

addpath(genpath('../../../Samak2.0'));

%% Test Spectra SAMAK and SSC 
A = InitKATRINDC('TD','Flat20');
samak = importdata('DC_Data/SpectraSAMAK/TBD_IntMatlab.txt');
ssc = importdata('DC_Data/SpectraSSC/TBD_IntSSC.txt');

samak(:,2) = samak(:,2) + 0.3;
ssc(:,2) = ssc(:,2) + 0.3;
samak(:,3) = sqrt(samak(:,2).*A.TimeSec.*A.qUfrac)./(A.TimeSec.*A.qUfrac);
ssc(:,3) = sqrt(ssc(:,2).*A.TimeSec.*A.qUfrac)./(A.TimeSec.*A.qUfrac);

comparespectra('SAMAKdata',samak,'SSCdata',ssc,'compType','INT',...
    'errorbarflag','ON','normFlag','ON')





