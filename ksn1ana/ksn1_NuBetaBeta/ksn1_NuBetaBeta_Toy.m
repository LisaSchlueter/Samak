%% load files
savedir = [getenv('Samak3.0'),'ksn1ana/ksn1_NuBetaBeta/results/'];
range = 40;
DataType = 'Real';
freePar = 'mNu E0 Norm Bkg';
pullFlag = 15;
switch pullFlag
    case 99
        pullStr = '';
    case 15
        pullStr = sprintf('_mNuSqPull1eV2');
end

file1 = sprintf('%sContour_Real_%s_40eV.mat',savedir,'E0NormBkg');
file2 = sprintf('%sContour_Real_%s_40eV.mat',savedir,'mNuE0NormBkg');
file3 = sprintf('%sContour_Real_%s_40eV_mNuSqPull1eV2.mat',savedir,'mNuE0NormBkg');

d1 = importdata(file1);
d2 = importdata(file2);
d3 = importdata(file3);

%%
mbbexp = @(phase,mbb3nu,sint4Sq,m4Sq) abs((1-sint4Sq).*mbb3nu+exp(1i*pi*phase).*sqrt(m4Sq).*sint4Sq);
Mode = 'NH';
if strcmp(Mode,'NH')
    % normal hierarchy: possible values in non-degenerate regime: 0-0.005 eV
    mbb3nu_min = 0;
    mbb3nu_max = 0.005;
elseif strcmp(Mode,'IH')
    % inverted hierarchy: possible values in non-degenerate regime: 0.01-0.05 eV
    mbb3nu_min    = 0.01;
    mbb3nu_max    = 0.05;
end

nGrid = 1e3;
sint4Sq = logspace(-3,log(0.75),nGrid);
m4Sq    = logspace(-1,3.2,nGrid);

%% extreme cases:
GetFigure;

contour(sint4Sq,m4Sq,mbbexp(1,mbb3nu_min,sint4Sq,m4Sq'),[0.2 0.2],'Color',rgb('DodgerBlue'));  % mbb3nu_min, +1
hold on; 
contour(sint4Sq,m4Sq,mbbexp(-1,mbb3nu_min,sint4Sq,m4Sq'),[0.2 0.2],'Color',rgb('Orange')); % mbb3nu_min, -1
contour(sint4Sq,m4Sq,mbbexp(0,mbb3nu_min,sint4Sq,m4Sq'),[0.2 0.2],'Color',rgb('IndianRed')); % mbb3nu_min, -0
contour(sint4Sq,m4Sq,mbbexp(1,mbb3nu_max,sint4Sq,m4Sq'),[0.2 0.2])
contour(sint4Sq,m4Sq,mbbexp(0,mbb3nu_max,sint4Sq,m4Sq'),[0.2 0.2])
contour(sint4Sq,m4Sq,mbbexp(-1,mbb3nu_max,sint4Sq,m4Sq'),[0.2 0.2])
%contour(sint4Sq,m4Sq,mbbexp(-1,mbb3nu_max,sint4Sq,m4Sq'),[0.2 0.2])
% contour(sint4Sq,m4Sq,mbbexp(0,mbb3nu_max,sint4Sq,m4Sq'),[0.2 0.2])
 set(gca,'YScale','log');
 set(gca,'XScale','log');
 hold on;
 pcontour = plot(d1.sinT4Sq,d1.mNu4Sq,'-','Color',rgb('Orange'),'LineWidth',2.5);
 hold off;
% hold off;
% 
