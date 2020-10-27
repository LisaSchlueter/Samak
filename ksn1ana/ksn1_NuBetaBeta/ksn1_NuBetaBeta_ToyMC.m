% settings
nGrid = 50;
nSamples = 5e3;
Mode = 'IH';
savedir = [getenv('SamakPath'),'ksn1ana/ksn1_NuBetaBeta/results/'];
mbbexplim = 0.165; % eV (KamLAND-Zen 90% C.L.)

% label
savename = sprintf('%sksn1_NuBetaBeta_ToyMC_%.3feVlim_%s_Grid%.0f_Samples%.0f_New.mat',savedir,mbbexplim,Mode,nGrid,nSamples);
if exist(savename,'file') 
    load(savename)
else

% grid
m4Sq    = repmat(logspace(-1,4,nGrid),nGrid,1);
sint4Sq = repmat(logspace(-3,log(0.74),nGrid),nGrid,1)';

phase    = 2.*rand(1,1,nSamples)-1;  % draw random phase for each grid point between +-1
if strcmp(Mode,'NH')
    % normal hierarchy: possible values in non-degenerate regime: 0-0.005 eV
    if contains(savename,'New')
        mbb3nu    = 0.026.*rand(1,1,nSamples); % using cosmology bound on m_light
    else
        mbb3nu    = 0.005.*rand(1,1,nSamples);
    end
    
elseif strcmp(Mode,'IH')
    % inverted hierarchy: possible values in non-degenerate regime: 0.01-0.05 eV
    if contains(savename,'New')
        mbb3nu    = (0.049-0.018).*rand(1,1,nSamples)+0.018; % using cosmology bound on m_light
    else
        mbb3nu    = 0.04.*rand(1,1,nSamples)+0.01;
    end
end

% mbb (4nu) formula
mbbexp = abs((1-sint4Sq).*mbb3nu+exp(1i.*pi.*phase).*sqrt(m4Sq).*sint4Sq);

sint2Sq_contour =repmat(logspace(-3,log(0.74),nGrid*50),nSamples,1)';%zeros(nGrid*10,nSamples);%cell(nSamples,1);
m4Sq_contour   = zeros(nGrid*50,nSamples);
Index = ones(nSamples,1);
for i=1:nSamples
    if all(all(mbbexp(:,:,i)<mbbexplim)) || all(all(mbbexp(:,:,i)>mbbexplim))
        Index(i) = 0;
        continue
    end
    [a,~] = contour(sint4Sq,m4Sq,squeeze(mbbexp(:,:,i)),[mbbexplim,mbbexplim]);
    hold on;
    try
        m4Sq_contour(:,i) = interp1(a(1,2:end),a(2,2:end),sint2Sq_contour(:,i),'lin','extrap');
    catch
        Index(i) = 0;
    end
end
close all
%% 
Index = logical(Index);
m4Sq_contour= m4Sq_contour(:,Index);
sint2Sq_contour = sint2Sq_contour(:,Index);
save(savename,'sint2Sq_contour','m4Sq_contour')
end

%% also save "light" version with less information
savename_light = strrep(savename,'.mat','_light.mat');
if ~exist(savename_light,'file') 
    sint2Sq     = sint2Sq_contour(:,1); 
    m4Sq_min     = min(m4Sq_contour');
    m4Sq_max     = max(m4Sq_contour');
    m4Sq_err     = 0.5*(m4Sq_max-m4Sq_min); % symmetric average err
    m4Sq      = m4Sq_max-m4Sq_err;        % middle
    save(savename_light,'sint2Sq','m4Sq','m4Sq_err','m4Sq_min','m4Sq_max');
else
    load(savename_light);
end
