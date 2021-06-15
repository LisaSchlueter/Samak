% Test of Wilk's theorem (coverage)
% debuging script
% make some sanity plots

%% load summary file
Hypothesis = 'H0';
SavePlt = 'ON';
NrandMC = 1e3;
switch Hypothesis
    case 'H0' 
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        yStr =  'H0';
    case 'H1'
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
         yStr =  'H1';
end
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_%.0fsamples.mat',savedir,NrandMC);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,NrandMC);
end

if exist(savefile,'file')
    load(savefile)
else
    return
end


%% chi2min vs chi2bf
nDelta = 8;
Colors = jet(nDelta);
pHandle = cell(nDelta,1);
legStr = cell(nDelta,1);

GetFigure;
x = linspace(0,max(chi2_null)+3,1e3);
pEq = plot(x,x,'-','LineWidth',2,'Color',rgb('Black'));
hold on;
for i=1:nDelta   
pHandle{i} = plot(chi2_bf(chi2_delta>(i-1) & chi2_delta<i),chi2_null(chi2_delta>(i-1) & chi2_delta<i),'.','MarkerSize',15,'Color',Colors(i,:));
legStr{i} = sprintf('%.0f < \\Delta\\chi^2 < %.0f',i-1,i);
end
PrettyFigureFormat('FontSize',22);
ylabel(sprintf('\\chi^2_{%s} (25 dof)',yStr));
xlabel(sprintf('\\chi^2_{min} (23 dof)'));

leg =legend([pEq,pHandle{:}],...
    sprintf('\\Delta\\chi^2 = 0'),...
    legStr{:},...
    'Location','southeast');
leg.Title.String = sprintf('\\Delta\\chi^2 = \\chi^2_{%s} - \\chi^2_{min}',yStr);
leg.Title.FontWeight = 'normal';
PrettyLegendFormat(leg);
xlim([0 max(x)])
ylim([0,max(chi2_bf)+3])
%% save
if strcmp(SavePlt,'ON')
 plotname = strrep(strrep(savefile,'results','plots'),'.mat','_Chi2minVsBfColor.png');
 print(gcf,plotname,'-dpng','-r450');
 fprintf('save plot to %s \n',plotname);
end
