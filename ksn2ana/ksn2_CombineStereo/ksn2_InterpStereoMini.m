% create mini file with interpolated STEREO data
% interpolated in log,log

% interpolate Stereo Chi^2 Map to match KATRIN grid
% sanity plots
Maxm4Sq    = 36^2;
Minm4Sq    = 0.1;
Type = 'Sensitivity';
RecomputeFlag = 'ON';

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/'];
savefile = sprintf('%sksn2_InterpStereoMini_Max%.0feV2_Min%.2gfeV2.mat',...
    savedir,Maxm4Sq,Minm4Sq);
if strcmp(Type,'Sensitivity')
    savefile = strrep(savefile,'.mat','_Sensitivity.mat');
end


if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    fprintf('file already exist %s \n',savefile);
else
    
    %% load original Stereo
    if strcmp(Type,'Sensitivity')
        StereoFile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/StereoMaps_Sensitivity.mat'];
    else
        StereoFile = [getenv('SamakPath'),'ksn2ana/ksn2_CombineStereo/results/StereoMaps.mat'];
    end

    ds = importdata(StereoFile);
    
    % KATRIN grid (interpolated)
    nInter      = 1e3; % always the same
    mNu4tmp     = logspace(log10(Minm4Sq),log10(Maxm4Sq),nInter);
    mNu4Sq      = repmat(mNu4tmp,nInter,1);
    sin2T4      = repmat(logspace(-3,log10(0.5),nInter),nInter,1)';
    sin2T4_Osci = 4*sin2T4.*(1-sin2T4);
    
    % find common  parameter space
    IntermNuIdx = mNu4Sq>=min(min(ds.mNu41Sq)) &  mNu4Sq<=max(max(ds.mNu41Sq));
    InterSinIdx = sin2T4_Osci>=min(min(ds.sin2TSq)) &  sin2T4_Osci<=max(max(ds.sin2TSq));
    InterIdx     = logical(IntermNuIdx.* InterSinIdx);
    [Startrow,Startcol]= find(InterIdx,1,'first');
    [Stoprow,Stopcol]= find(InterIdx,1,'last');
    
    %interpolte STEREO
    mNu4Sq_cut = mNu4Sq(Startrow:Stoprow,Startcol:Stopcol);%repmat(mNu4Sq_tmp,numel(sin2T4_tmp),1); % part of KATRIN that is interpolated
    sin2T4_cut = sin2T4(Startrow:Stoprow,Startcol:Stopcol);% repmat(sin2T4_tmp,1,numel(mNu4Sq_tmp));
    sin2T4_Osci_cut = sin2T4_Osci(Startrow:Stoprow,Startcol:Stopcol);%
    chi2Stereo_cut     = interp2(ds.mNu41Sq,ds.sin2TSq,ds.chi2,mNu4Sq_cut,sin2T4_Osci_cut,'spline');
    chi2Stereocrit_cut = interp2(ds.mNu41Sq,ds.sin2TSq,ds.chi2crit,mNu4Sq_cut,sin2T4_Osci_cut,'spline');
    chi2Stereo_min = min(min(chi2Stereo_cut)); % not zero, because best-fit not on grid point
    
    % intert into KATRIN par space
    DeltaChi2Wilks = GetDeltaChi2(95,2);
    if strcmp(Type,'Sensitivity')
        %sensitivity delta-chi^2 map
        chi2Stereo = chi2Stereo_min.*ones(nInter,nInter);
    else
        chi2Stereo = zeros(nInter,nInter);
    end
    
    chi2Stereo(InterIdx) = chi2Stereo_cut;
    chi2Stereocrit = DeltaChi2Wilks.*ones(nInter,nInter);
    chi2Stereocrit(InterIdx) = chi2Stereocrit_cut;
    %%
   
    
    save(savefile,'mNu4Sq_cut','sin2T4_cut','sin2T4_Osci_cut','chi2Stereo_cut','chi2Stereocrit_cut',...
        'mNu4Sq','sin2T4','sin2T4_Osci','chi2Stereo','chi2Stereocrit',...
        'InterIdx','Startrow','Startcol','Stoprow','Stopcol');
    fprintf('save to %s \n',savefile)
end


