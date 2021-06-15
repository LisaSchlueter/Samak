function d = ksn2_GetRandTritiumSpectrum(Hypothesis)
% short cut function to load random tritium spectrum file
%% load random tritium spectra
range = 40;
 chi2 = 'chi2CMShape';
 
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if strcmp(Hypothesis,'H0')
         Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;       
elseif strcmp(Hypothesis,'H1')
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
end

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_RandomSpectra_%s_%.0feV_NullHypothesis_%.0fsamples.mat',savedir,chi2,range,5000);
else
    savefile = sprintf('%sksn2_WilksTheorem_RandomSpectra_%s_%.0feV_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',...
        savedir,chi2,range,Twin_mNu4Sq,Twin_sin2T4,5000);
end

if exist(savefile,'file')
    d = importdata(savefile);
else
    fprintf('file %s not found \n',savefile)
    return
end
end