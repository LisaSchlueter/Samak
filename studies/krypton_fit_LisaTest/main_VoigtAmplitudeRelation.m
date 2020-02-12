clear all;
addpath(genpath('../../../Samak2.0'));

temperature = 1:10:351;
amplitude = [];
for i= 1:36

myKrOpt= {'TD', 'KrL3_32_Satellites','FPD_Segmentation','OFF' , ...
    'HVRipples','OFF', 'MultiPeaksFlag','OFF',...
    'DopplerEffectFlag','ON', 'ConvFlag','DopplerEffect', ...
    'T', temperature(i)};
myObj = InitKrKATRIN_LisaFit(myKrOpt{:});
myObj.ComputeKrDS();
amplitude_voigt(i) = max(myObj.KrDS);
end 

myObjNoDoppler = InitKrKATRIN_LisaFit('TD', 'KrL3_32_Satellites','FPD_Segmentation','OFF' , ...
    'HVRipples','OFF', 'MultiPeaksFlag','OFF',...
    'DopplerEffectFlag','OFF', 'ConvFlag','OFF');
myObjNoDoppler.ComputeKrDS();
amplitude_nodopper = max(myObjNoDoppler.KrDS);
figure(99);
%% 
plot(temperature, amplitude_voigt,'-');
