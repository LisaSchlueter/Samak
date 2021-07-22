% calculate 0NuBetaBeta contour numerically
% most stringent contour always for phase = 0 and maximal m3nubetabeta
% --> term with m4 gets maximal -> m4 need to be even smaller to compensate

% least stringend contour phase = +-1 and  maximal m3nubetabeta 

%% m_betabeta neutrino mass in 3+1 neutrino framework
mbb4nu = @(phase,mbb3nu,sint4Sq,m4Sq) abs((1-sint4Sq).*mbb3nu+exp(1i*pi*phase).*sqrt(m4Sq).*sint4Sq);

%% get lines
nGrid = 1e3;
sin2T4 = repmat(logspace(-3,log10(1),nGrid)',[1,nGrid]);% mixing
m4Sq   = repmat(logspace(-1,log10(40^2),nGrid),[nGrid,1]);

mbb3nu_min = [0,0.018]; 
mbb3nu_max = [0.005,0.05];%[0.026,0.049]';%NH, IH
Phase = [0,1];

mbb4nu_ExpLimit =  0.165;

%% Get all the contours: NH
Idx = 1;
mbb4nu_NH_0 = mbb4nu(0,mbb3nu_max(Idx),sin2T4,m4Sq);
mbb4nu_NH_1max = mbb4nu(1,mbb3nu_max(Idx),sin2T4,m4Sq);
mbb4nu_NH_1min = mbb4nu(1,mbb3nu_min(Idx),sin2T4,m4Sq);
mbb4nu_NH_m1max = mbb4nu(-1,mbb3nu_max(Idx),sin2T4,m4Sq);
mbb4nu_NH_m1min = mbb4nu(-1,mbb3nu_min(Idx),sin2T4,m4Sq);
% Get all the contours
[Contour_NH_0,~] = contour(sin2T4,m4Sq,mbb4nu_NH_0,mbb4nu_ExpLimit.*ones(2,1));
[Contour_NH_1max,~] = contour(sin2T4,m4Sq,mbb4nu_NH_1max,mbb4nu_ExpLimit.*ones(2,1));
[Contour_NH_1min,~] = contour(sin2T4,m4Sq,mbb4nu_NH_1min,mbb4nu_ExpLimit.*ones(2,1));
[Contour_NH_m1max,~] = contour(sin2T4,m4Sq,mbb4nu_NH_m1max,mbb4nu_ExpLimit.*ones(2,1));
[Contour_NH_m1min,~] = contour(sin2T4,m4Sq,mbb4nu_NH_m1min,mbb4nu_ExpLimit.*ones(2,1));
close all;

%% Get all the contours: IH
Idx = 2;
mbb4nu_IH_0 = mbb4nu(0,mbb3nu_max(Idx),sin2T4,m4Sq);
mbb4nu_IH_1max = mbb4nu(1,mbb3nu_max(Idx),sin2T4,m4Sq);
mbb4nu_IH_1min = mbb4nu(1,mbb3nu_min(Idx),sin2T4,m4Sq);
mbb4nu_IH_m1max = mbb4nu(-1,mbb3nu_max(Idx),sin2T4,m4Sq);
mbb4nu_IH_m1min = mbb4nu(-1,mbb3nu_min(Idx),sin2T4,m4Sq);
% Get all the contours
[Contour_IH_0,~] = contour(sin2T4,m4Sq,mbb4nu_IH_0,mbb4nu_ExpLimit.*ones(2,1));
[Contour_IH_1max,~] = contour(sin2T4,m4Sq,mbb4nu_IH_1max,mbb4nu_ExpLimit.*ones(2,1));
[Contour_IH_1min,~] = contour(sin2T4,m4Sq,mbb4nu_IH_1min,mbb4nu_ExpLimit.*ones(2,1));
[Contour_IH_m1max,~] = contour(sin2T4,m4Sq,mbb4nu_IH_m1max,mbb4nu_ExpLimit.*ones(2,1));
[Contour_IH_m1min,~] = contour(sin2T4,m4Sq,mbb4nu_IH_m1min,mbb4nu_ExpLimit.*ones(2,1));
close all;

%%
GetFigure;
p0IH = plot(Contour_IH_0(1,2:end),Contour_IH_0(2,2:end),'k-','LineWidth',2);
hold on;
p1maxIH = plot(Contour_IH_1max(1,2:end),Contour_IH_1max(2,2:end),'-','Color',rgb('Orange'),'LineWidth',2);
p1minIH = plot(Contour_IH_1min(1,2:end),Contour_IH_1min(2,2:end),'-','Color',rgb('IndianRed'),'LineWidth',2);


p0NH = plot(Contour_NH_0(1,2:end),Contour_NH_0(2,2:end),':','Color',rgb('Silver'),'LineWidth',2);
p1maxNH = plot(Contour_NH_1max(1,2:end),Contour_NH_1max(2,2:end),':','Color',rgb('DodgerBlue'),'LineWidth',2.5);
p1minNH = plot(Contour_NH_1min(1,2:end),Contour_NH_1min(2,2:end),':','Color',rgb('Navy'),'LineWidth',2.5);

set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([6e-03,1]);
ylim([0.1,1600]);
PrettyFigureFormat;
xlabel(sprintf('|{\\it U_{e4}}|^2'));
ylabel(sprintf('{\\it m}_4^2 (eV^2)'));

leg = legend([p0IH,p1maxIH,p1minIH,p0NH,p1maxNH,p1minNH],...
    'IH Phase 0','IH Phase 1 , max mbb3nu','IH Phase 1 , min mbb3nu',...
    'NH Phase 0','NH Phase 1 , max mbb3nu','NH Phase 1 , min mbb3nu',...
        'box','off');
    leg.NumColumns = 2;
