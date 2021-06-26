% read original Stereo sensitivity
% transform into mat file for later use

matfile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/StereoMaps_Sensitivity.mat'];

if exist(matfile,'file')
    load(matfile)
else 
    Sfile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/MapDeltaChisquareSensitivityValues.txt'];
    Scritfile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/Stereo_Chi2CritMap.csv'];
    
    d = importdata(Sfile);
    sin2TSq = d.data(:,1); sin2TSq = reshape(sin2TSq,100,90);
    mNu41Sq = d.data(:,2); mNu41Sq = reshape(mNu41Sq,100,90);
    chi2    = d.data(:,3); chi2    = reshape(chi2,100,90);
    
    datacrit = csvread(Scritfile,13,0,[13 0 9012 2]);
    chi2crit    = datacrit(:,3); chi2crit    = reshape(chi2crit,100,90);
    
    % make sure same m4,sin are the same
    if any(datacrit(:,1)~=d.data(:,1)) ||  any(datacrit(:,2)~=d.data(:,2))
        fprintf(2,'chi2 and chi2-crit dont have the same binning \n')
        return
    end
    save(matfile,'sin2TSq','mNu41Sq','chi2','chi2crit');
    fprintf('save file to %s \n',matfile)
end
