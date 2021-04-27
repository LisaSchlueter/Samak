% plot closed contours
Hypothesis = 'H0';
switch Hypothesis
    case 'H0'
        randMC = 1:1000;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
    case 'H1'
        randMC = [1:352,384:840];
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2Stat';
end
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_%.0fsamples.mat',savedir,numel(randMC));
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,numel(randMC));
end

if exist(savefile,'file') 
    fprintf('savefile already created \n');
    d = importdata(savefile);
else
    return
end

nClosed = sum(d.ClosedLog95);
IdxClosed = find(d.ClosedLog95);
PltColors = jet(nClosed);
 
% GetFigure;
% for i=1:nClosed
% %plot(d.sin2T4_contour{IdxClosed(i)},d.mNu4Sq_contour{IdxClosed(i)},'.','Color',PltColors(i,:));
% plot(d.sin2T4_bf(IdxClosed(i)),d.mNu4Sq_bf(IdxClosed(i)),'x',...
%     'Color',PltColors(i,:),'MarkerSize',10,'LineWidth',1.5);
% 
% hold on;
% 
% end
% set(gca,'YScale','log');
% set(gca,'XScale','log');


%%
GetFigure
PltColors2 = jet(1e3);
pSens = plot(d.sin2T4_contour_Asimov,d.mNu4Sq_contour_Asimov,'-k','LineWidth',2);
hold on;
x = linspace(1e-03,0.5,1e2);
y = linspace(0.1,40^2,1e2);
pref = plot(x,0.1.*ones(100,1),'-','Color',rgb('Silver'),'LineWidth',1.5);
plot(x,40^2.*ones(100,1),'-','Color',rgb('Silver'),'LineWidth',1.5);
plot(1e-03.*ones(100,1),y,'-','Color',rgb('Silver'),'LineWidth',1.5);
plot(0.5.*ones(100,1),y,'-','Color',rgb('Silver'),'LineWidth',1.5);

% look if best fit is in our outside
InOutIdx = zeros(1e3,1);

for i=1:1e3
     p1 = plot(d.sin2T4_bf(i),d.mNu4Sq_bf(i),'.',...
        'Color',rgb('Orange'),'MarkerSize',10,'LineWidth',1.5);

   if i~=IdxClosed
        p1.Color = rgb('SkyBlue');
        p1.LineWidth = 1;
        pHandleOpen = p1;
   else 
       pHandleReg = p1;
     %  p1.Marker = '.'; %p1.MarkerSize = 4;
       p1.MarkerFaceColor = rgb('Orange');
   end  
   
   %look if best fit is in our outside
   sinTmp = interp1(d.mNu4Sq_contour_Asimov,d.sin2T4_contour_Asimov,d.mNu4Sq_bf(i),'spline');
   InOutIdx(i) = sinTmp<=d.sin2T4_bf(i); % 1 means inside (significant=
   if InOutIdx(i)
       p1.Color = rgb('IndianRed');
       pHandleIn = p1;
   end
end

set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);
leg = legend([pSens,pref,pHandleOpen,pHandleReg,pHandleIn],...
    'Sensitivity','Grid boundaries',sprintf('Rand MC: \\Delta\\chi^2 < 5.99 (%.1f%%)',100.*(numel(d.sin2T4_bf)-numel(IdxClosed))./numel(d.sin2T4_bf)),...
    sprintf('Rand MC: \\Delta\\chi^2 \\geq 5.99 (%.1f%%)',100*numel(IdxClosed)./numel(d.sin2T4_bf)),...
    sprintf('Rand MC: Best fit inside Sensitivity curve (%.1f%%)',100.*sum(InOutIdx)./numel(d.sin2T4_bf)),...
    'Location','southwest');
PrettyLegendFormat(leg);
xlim([8e-04 0.65]);
ylim([0.07 50^2]);

plotnameContourBf = strrep(strrep(savefile,'results','plots'),'.mat','_BestFitPosition.png');
print(gcf,plotnameContourBf,'-dpng','-r450');
fprintf('save plot to %s \n',plotnameContourBf);
%%
% Idx = d.sin2T4_bf~=0.5;
% ClosedContour = chi2_delta(Idx)>=5.99;
% sum(ClosedContour)./numel(ClosedContour)*100;
GetFigure;
h1 = histogram(d.chi2_null);
hold on;
h2 = histogram(d.chi2_null(d.ClosedLog95),'BinWidth',h1.BinWidth);
xlabel('chi2 null')

GetFigure;
h1 = histogram(d.chi2_bf);
hold on;
h2 = histogram(d.chi2_bf(d.ClosedLog95),'BinWidth',h1.BinWidth);
xlabel('chi2 bf')
