% Read KNM1 Data
[Real , Twin] = Knm1RealTwin_Create('Real','ON','Twin','ON');

%% Tritium Activity
Real.PlotSCdata_TritiumActivity;

%% Subrun Activity Fluctuation 
myMainTitle = sprintf('KNM1 subrun activity fluctuations - %.0f runs - %.0f subruns',numel(Real.RunList),numel(Real.RunList)*numel(Real.RunData.qU));
maintitle   = myMainTitle;
savefile    = sprintf('plots/KNM1_subrunactivityfluctuations_1.png');
fig1      = figure('Name','KNM1 subrun activity fluctuations - 274 runs - 10960 subruns','NumberTitle','off','rend','painters','pos',[10 10 1200 1200]);
a=annotation('textbox', [0 0.91 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=22;a.FontWeight='bold';
StackedSubRunWiseSTD=zeros(1,numel(Real.RunData.qU));
for qUindex = 1:numel(Real.RunData.qU)
SubRunActivity = ...
    (Real.SingleRunData.WGTS_MolFrac_TT_SubRun(qUindex,:)        + ...
    0.5 * Real.SingleRunData.WGTS_MolFrac_HT_SubRun(qUindex,:)   + ...
    0.5 * Real.SingleRunData.WGTS_MolFrac_DT_SubRun(qUindex,:)) .* ...
    Real.SingleRunData.WGTS_CD_MolPerCm2_SubRun(qUindex,:)      .* ...
    Real.ModelObj.ISXsection*1e4;
SubRunActivity=SubRunActivity(SubRunActivity>0);
SubRunActivity=SubRunActivity./mean(SubRunActivity);

% Error on Mean Activity for stacked spectrum
% http://davidmlane.com/hyperstat/A103735.html
StackedSubRunWiseSTD(qUindex)=std(SubRunActivity)/sqrt(numel(Real.RunList))*100;

fprintf('qU=%.0f - mean = %.5f sigma = %.5f \n',qUindex,mean(SubRunActivity),std(SubRunActivity));
subplot(numel(Real.RunData.qU)/4,numel(Real.RunData.qU)/10,qUindex)
histogram(SubRunActivity(SubRunActivity>0),'Normalization','Probability','LineWidth',2);
title(sprintf('subrun %.0f',qUindex));
  PrettyFigureFormat
  set(gca,'FontSize',12);
hold on
end
hold off
export_fig(gcf,savefile,'-q101','-m3');

%%
%% Stacked SubRun wise Activity Fluctuation 
myMainTitle = sprintf('KNM1 subrun activity fluctuations for %.0f runs stacked',numel(Real.RunList));
maintitle   = myMainTitle;
savefile    = sprintf('plots/KNM1_subrunactivityfluctuations_2.png');
fig1      = figure('Name','KNM1 subrun activity fluctuations: stacked','NumberTitle','off','rend','painters','pos',[10 10 1400 600]);
a=annotation('textbox', [0 0.91 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=22;a.FontWeight='bold';
subplot(1,4,[1 3])
plot(1:numel(Real.RunData.qU),StackedSubRunWiseSTD,'ks-.','MarkerSize',8,'MarkerFaceColor',rgb('SteelBlue'),'LineWidth',1);
ylabel('TASR Relative Uncertainty(%)');
xlabel('Subrun Index');
PrettyFigureFormat; set(gca,'FontSize',20);
subplot(1,4,[4])
histogram(StackedSubRunWiseSTD,4,'Normalization','Probability','LineWidth',2);
xlabel('TASR Relative Uncertainty(%)');
ylabel('stacked spectrum subrun frequency');
PrettyFigureFormat; set(gca,'FontSize',20);
 export_fig(gcf,savefile,'-q101','-m3');
 
