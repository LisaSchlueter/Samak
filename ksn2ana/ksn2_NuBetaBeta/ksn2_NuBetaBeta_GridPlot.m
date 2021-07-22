% calculate 0NuBetaBeta contour
% using grid 

% m_betabeta neutrino mass in 3+1 neutrino framework
mbb4nu = @(phase,mbb3nu,sint4Sq,m4Sq) abs((1-sint4Sq).*mbb3nu+exp(1i*pi*phase).*sqrt(m4Sq).*sint4Sq);
mbb3nu_min = [0,0.018]; 
mbb3nu_max = [0.026,0.049]';%NH, IH
Phase = [0,1];

nGrid = 1e3;
sin2T4 = repmat(logspace(-3,log10(1),nGrid)',[1,nGrid]);% mixing
m4Sq   = repmat(logspace(-2,log10(40^2),nGrid),[nGrid,1]);

mbb4nu_IH0  = mbb4nu(0,mbb3nu_max(2),sin2T4,m4Sq);
mbb4nu_IH1  = mbb4nu(1,mbb3nu_max(2),sin2T4,m4Sq);
mbb4nu_IHm1 = mbb4nu(-1,mbb3nu_max(2),sin2T4,m4Sq);
mbb4nu_ExpLimit =  0.165;


close all;
GetFigure;
[M,pc] = contour(sin2T4,m4Sq,mbb4nu_IH0,[mbb4nu_ExpLimit mbb4nu_ExpLimit]);
pc.delete;
mbb4nu_IH0(mbb4nu_IH0>mbb4nu_ExpLimit) = NaN;
surf(sin2T4,m4Sq,mbb4nu_IH0,'EdgeColor','none','FaceColor','interp');
hold on;
plot3(M(1,2:end),M(2,2:end),ones(size(M,2)-1,1),'k-','LineWidth',2);
view(2)
grid off;
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([6e-03,1]);
ylim([0.1,1600]);
c = colorbar; 
c.Label.String = sprintf('{\\it m}_{\\beta\\beta}^{4\\nu} (eV)');
xlabel(sprintf('|{\\it U_{e4}}|^2'));
ylabel(sprintf('{\\it m}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);
