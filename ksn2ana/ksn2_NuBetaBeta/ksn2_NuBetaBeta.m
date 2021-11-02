% calculate 0NuBetaBeta contour numerically
% most stringent contour always for phase = 0 and maximal m3nubetabeta
% --> term with m4 gets maximal -> m4 need to be even smaller to compensate
% least stringend contour phase = +-1 and  maximal m3nubetabeta
mbb3nu_max = [0.005,0.05];%[0.026,0.049]';%NH, IH
mbb4nu_ExpLimit =  0.16;%0.165;

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_NuBetaBeta/results/'];
savefile = sprintf('%sksn2_NuBetaBeta_mbb3numax_NH%.3feV_IH%.3feV_mbbExpLimit%.3feV.mat',...
    savedir,mbb3nu_max(1),mbb3nu_max(2),mbb4nu_ExpLimit);
if exist(savefile,'file') 
    load(savefile)
else
    %% m_betabeta neutrino mass in 3+1 neutrino framework
    mbb4nu = @(phase,mbb3nu,sint4Sq,m4Sq) abs((1-sint4Sq).*mbb3nu+exp(1i*pi*phase).*sqrt(m4Sq).*sint4Sq);
    
    nGrid = 1e3;
    sin2T4 = repmat(logspace(-3,log10(2),nGrid)',[1,nGrid]);% mixing
    m4Sq   = repmat(logspace(-1.58,log10(8e3),nGrid),[nGrid,1]);
    
    %% extreme cases for IH:
    mbb4nu_IH_min = mbb4nu(0,mbb3nu_max(2),sin2T4,m4Sq);
    mbb4nu_IH_max = mbb4nu(1,mbb3nu_max(2),sin2T4,m4Sq);
    
    % Get all the contours
    [Contour_IH_min,~] = contour(sin2T4,m4Sq,mbb4nu_IH_min,mbb4nu_ExpLimit.*ones(2,1));
    [Contour_IH_max,~] = contour(sin2T4,m4Sq,mbb4nu_IH_max,mbb4nu_ExpLimit.*ones(2,1));
    close all;
    
    % interpolate
    sin2T4_IH = logspace(log10(max([min(Contour_IH_min(1,2:end)),min(Contour_IH_max(1,2:end))])),...
        log10(min([max(Contour_IH_min(1,2:end)),max(Contour_IH_max(1,2:end))])),...
        1e3);
    
    m4Sq_IH_min = interp1(Contour_IH_min(1,2:end),Contour_IH_min(2,2:end),sin2T4_IH,'spline');
    m4Sq_IH_max = interp1(Contour_IH_max(1,2:end),Contour_IH_max(2,2:end),sin2T4_IH,'spline');
    m4Sq_IH_mean = mean([m4Sq_IH_max;m4Sq_IH_min]);
    
    %% extreme cases for NH:
    mbb4nu_NH_min = mbb4nu(0,mbb3nu_max(1),sin2T4,m4Sq);
    mbb4nu_NH_max = mbb4nu(1,mbb3nu_max(1),sin2T4,m4Sq);
    
    % Get all the contours
    [Contour_NH_min,~] = contour(sin2T4,m4Sq,mbb4nu_NH_min,mbb4nu_ExpLimit.*ones(2,1));
    [Contour_NH_max,~] = contour(sin2T4,m4Sq,mbb4nu_NH_max,mbb4nu_ExpLimit.*ones(2,1));
    close all;
    
    % interpolate
    sin2T4_NH = logspace(log10(max([min(Contour_NH_min(1,2:end)),min(Contour_NH_max(1,2:end))])),...
        log10(min([max(Contour_NH_min(1,2:end)),max(Contour_NH_max(1,2:end))])),...
        1e3);
    m4Sq_NH_min = interp1(Contour_NH_min(1,2:end),Contour_NH_min(2,2:end),sin2T4_NH,'spline');
    m4Sq_NH_max = interp1(Contour_NH_max(1,2:end),Contour_NH_max(2,2:end),sin2T4_NH,'spline');
    m4Sq_NH_mean = mean([m4Sq_NH_max;m4Sq_NH_min]);
    
    MakeDir(savedir);
    save(savefile,'sin2T4_IH','m4Sq_IH_min','m4Sq_IH_max','m4Sq_IH_mean',...
        'sin2T4_NH','m4Sq_NH_min','m4Sq_NH_max','m4Sq_NH_mean',...
        'mbb3nu_max','mbb4nu_ExpLimit');
    fprintf('save file to %s \n',savefile);
end

%% plot result
xIH = 4.*sin2T4_IH.*(1-sin2T4_IH);
xNH = 4.*sin2T4_NH.*(1-sin2T4_NH);
GetFigure

[l,aIH]= boundedline(xIH,m4Sq_IH_mean,abs(m4Sq_IH_mean-[m4Sq_IH_max;m4Sq_IH_min])');
l.delete; aIH.FaceColor = rgb('LightGray');
hold on;
% plot(xIH,m4Sq_IH_min,'k-','LineWidth',2)
% plot(xIH,m4Sq_IH_max,'k-','LineWidth',2)

[l,aNH]= boundedline(xNH,m4Sq_NH_mean,abs(m4Sq_NH_mean-[m4Sq_NH_max;m4Sq_NH_min])');
l.delete;
aNH.FaceColor = rgb('IndianRed');
% plot(xNH,m4Sq_NH_min,'-','LineWidth',2,'Color',rgb('FireBrick'));
% plot(xNH,m4Sq_NH_max,'-','LineWidth',2,'Color',rgb('FireBrick'));

set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([6e-03,1]);
ylim([0.1,1600]);
PrettyFigureFormat;
xlabel(sprintf('sin^2(2\\theta_{ee})'));%|{\\itU_{e4}}|^2'));
ylabel(sprintf('{\\it m}_4^2 (eV^2)'));
leg = legend([aNH,aIH],'Normal Ordering','Inverted Ordering','box','off');
