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
% chi2_null = ReCalc_chi2Null_i;
GetFigure;
x = linspace(0,max(chi2_null)+3,1e3);
pEq = plot(x,x,'-','LineWidth',2,'Color',rgb('Black'));
hold on;
if strcmp(Hypothesis,'H1')
    pSig = plot(chi2_bf(ClosedLog95),chi2_null(ClosedLog95),'.','MarkerSize',15,'Color',rgb('Orange'));
else
   pSig = plot(chi2_bf(ClosedLog95),chi2_null(ClosedLog95),'.','MarkerSize',15,'Color',rgb('Orange'));
   InOutIdx = logical(InOutIdx);
   pSig2 = plot(chi2_bf(InOutIdx),chi2_null(InOutIdx),'.','MarkerSize',15,'Color',rgb('IndianRed'));  
end
pNonSig = plot(chi2_bf(~ClosedLog95),chi2_null(~ClosedLog95),'.','MarkerSize',15,'Color',rgb('SkyBlue'));
PrettyFigureFormat('FontSize',22);
ylabel(sprintf('\\chi^2_{%s} (25 dof)',yStr));
xlabel(sprintf('\\chi^2_{min} (23 dof)'));
%%
if strcmp(Hypothesis,'H1')
leg =legend([pEq,pNonSig,pSig],...
    sprintf('\\Delta\\chi^2 = \\chi^2_{%s} - \\chi^2_{min} = 0',yStr),...
    sprintf('\\Delta\\chi^{2 } < 5.99'),...
    sprintf('\\Delta\\chi^2 \\geq 5.99'),'Location','northwest');
else
    leg =legend([pEq,pNonSig,pSig,pSig2],...
    sprintf('\\Delta\\chi^2 = \\chi^2_{%s} - \\chi^2_{min} = 0',yStr),...
    sprintf('\\Delta\\chi^{2 } < 5.99 (%.1f%%)',1e2*sum(~ClosedLog95)./NrandMC),...
     sprintf('\\Delta\\chi^2 \\geq 5.99 (%.1f%%)',1e2*sum(ClosedLog95)./NrandMC),...
     sprintf('\\Delta\\chi^2 \\geq 5.99 & inside sensitivity curve (%.1f%%)',1e2.*sum(InOutIdx)./NrandMC),...
     'Location','northwest');
end
PrettyLegendFormat(leg);
xlim([0 max(x)])

%% save
if strcmp(SavePlt,'ON')
 plotname = strrep(strrep(savefile,'results','plots'),'.mat','_Chi2minVsBf.png');
 print(gcf,plotname,'-dpng','-r450');
 fprintf('save plot to %s \n',plotname);
end
