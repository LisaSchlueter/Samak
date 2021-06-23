% read original Stereo data
% transform into mat file for later use

matfile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/StereoMaps.mat'];
if exist(matfile,'file') &&1==2
    load(matfile)
else
    Sfile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/Stereo_Chi2Map.csv'];
    Scritfile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/Stereo_Chi2CritMap.csv'];
    
    data = csvread(Sfile,8,0);
    sin2TSq = data(:,1); sin2TSq = reshape(sin2TSq,100,90);
    mNu41Sq = data(:,2); mNu41Sq = reshape(mNu41Sq,100,90);
    chi2    = data(:,3); chi2    = reshape(chi2,100,90);
    
    datacrit = csvread(Scritfile,13,0,[13 0 9012 2]);
    chi2crit    = datacrit(:,3); chi2crit    = reshape(chi2crit,100,90);
    
    % make sure same m4,sin are the same
    if any(datacrit(:,1)~=data(:,1)) ||  any(datacrit(:,2)~=data(:,2))
        return
    end
    save(matfile,'sin2TSq','mNu41Sq','chi2','chi2crit');
end
