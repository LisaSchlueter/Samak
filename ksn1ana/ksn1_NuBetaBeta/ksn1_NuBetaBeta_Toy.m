%% load files
savedir = [getenv('SamakPath'),'ksn1ana/ksn1_NuBetaBeta/results/'];
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

%% convert sint^2 to sin2t^2
sint4Sq_Osci       = Convert2Osci(sint4Sq);
sint4SqKATRIN_Osci = Convert2Osci(d1.sinT4Sq);
%% special cases
Phase = [-1,-0.5,0,0.5,1];
m4Sq_contourMin = zeros(numel(m4Sq),numel(Phase));
m4Sq_contourMax = zeros(numel(m4Sq),numel(Phase));
for i=1:numel(Phase)
    [a,] = contour(sint4Sq,m4Sq,mbbexp(Phase(i),mbb3nu_min,sint4Sq,m4Sq'),[0.165 0.165],'Color',rgb('DodgerBlue'));  % mbb3nu_min, +1
    m4Sq_contourMin(:,i) = interp1(a(1,2:end),a(2,2:end),sint4Sq,'lin','extrap');
    
    close;
    
    [a,] = contour(sint4Sq,m4Sq,mbbexp(Phase(i),mbb3nu_max,sint4Sq,m4Sq'),[0.165 0.165],'Color',rgb('DodgerBlue'));  % mbb3nu_min, +1
    m4Sq_contourMax(:,i) = interp1(a(1,2:end),a(2,2:end),sint4Sq,'lin','extrap');
    
    close;
end
%% plot some special cases:
MyColors = hsv(numel(Phase));
f1 = figure('Units','normalized','Position',[0.1,0.1,0.8,0.4]);
subplot(1,2,1);
p = cell(numel(Phase),1);
legStr = cell(numel(Phase),1);
for i=1:numel(Phase)
    % p{i}= plot(sint4Sq,m4Sq_contourMin(:,i),'-','Color',rgb('DodgerBlue'));
    
    if Phase(i)<0
        LineStyle = '--';
    elseif Phase(i)>0
        LineStyle = ':';
    elseif Phase(i)==0
        LineStyle = '-';
    end
    
    p{i}= plot(sint4Sq_Osci,m4Sq_contourMin(:,i),'LineStyle',LineStyle,'Color',MyColors(i,:),'LineWidth',1);
    hold on;
    legStr{i} = sprintf('m_{\\beta\\beta}^{3\\nu} = %.3g eV, \\phi = %.1g \\pi',mbb3nu_min,Phase(i));
end

pcontour = plot(sint4SqKATRIN_Osci,d1.mNu4Sq,'-','Color',rgb('Black'),'LineWidth',1.5);
set(gca,'XScale','log');
set(gca,'YScale','log');
PrettyFigureFormat;

leg = legend([pcontour,p{:}],{'KATRIN 95% C.L.',legStr{:}},'Location','southwest');

legend box off;
xlabel(sprintf('sin^2(2{\\it\\theta}_{ee})'));ylabel(sprintf('\\Delta m_{14}^2 (eV^2)'));
ylabel(sprintf('\\Delta m_{14}^2 (eV^2)'));

subplot(1,2,2);
p = cell(numel(Phase),1);
for i=1:numel(Phase)
   if Phase(i)<0
        LineStyle = '--';
    elseif Phase(i)>0
        LineStyle = ':';
    elseif Phase(i)==0
        LineStyle = '-';
    end
    p{i}= plot(sint4Sq_Osci,m4Sq_contourMax(:,i),'LineStyle',LineStyle,'Color',MyColors(i,:),'LineWidth',1.5);
    hold on;
    legStr{i} = sprintf('m_{\\beta\\beta}^{3\\nu} = %.3g eV, \\phi = %.1g \\pi',mbb3nu_max,Phase(i));
end

pcontour = plot(sint4SqKATRIN_Osci,d1.mNu4Sq,'-','Color',rgb('Black'),'LineWidth',1.5);
set(gca,'XScale','log');
set(gca,'YScale','log');
PrettyFigureFormat;

leg = legend([pcontour,p{:}],{'KATRIN 95% C.L.',legStr{:}},'Location','southwest');
legend box off;

t = sgtitle(sprintf('%s, m_{\\beta\\beta}^{3\\nu} \\in [%.3g,%.3g] eV , \\phi \\in [-\\pi,\\pi] ',Mode,mbb3nu_min,mbb3nu_max));
t.FontSize = 18;

xlabel(sprintf('sin^2(2{\\it\\theta}_{ee})'));
ylabel(sprintf('\\Delta m_{14}^2 (eV^2)'));
%%
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = sprintf('%sksn1_NuBetaBeta_ToyMC_%s',plotdir,Mode);
print(plotname,'-dpng','-r300');

%%
function out = Convert2Osci(sinT4Sq)
out = 4*sinT4Sq.*(1-sinT4Sq);
end