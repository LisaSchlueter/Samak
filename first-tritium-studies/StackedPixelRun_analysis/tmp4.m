myindex =  -[1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];
myindex=myindex(1:18);

ab_e0_stat    = [18573.553      18573.5808      18573.6267      18573.6331      18573.5939       18573.578      18573.5655       18573.555      18573.7422      18573.8384      18573.8853      18573.8563      18574.2124      18574.1352      18574.3116      18574.2896      18574.4059      18574.0543];
ab_e0err_stat = [0.02721    0.028956    0.030925     0.03328    0.036349    0.041032    0.050318    0.061547    0.085396    0.098246     0.11842     0.15036      0.1969     0.22922     0.27391     0.33368     0.41173     0.47709];
 
cw_e0_stat    = [18572.9342      18573.1274      18573.3297      18573.4803      18573.5721      18573.6754      18573.7675      18573.7866      18573.9708      18574.0554      18574.0834      18574.0259      18574.3471      18574.2504      18574.4056      18574.3577      18574.4447      18574.0624];
cw_e0err_stat = [0.02748    0.029191    0.031114    0.033419    0.036438    0.041077     0.05033    0.061569    0.085486    0.098377     0.11861     0.15066     0.19743      0.2299     0.27482     0.33489     0.41326     0.47835];

mode='Data';
fixPar = '1 5 6';
FixParLabel = fixPar(~isspace(fixPar));
RunList = 'StackCD100all';

%% Endpoint Verus Range
myMainTitle=[sprintf('KATRIN First Tritium Samak Fit - %s - qU Scan - FixPar = %s',mode,fixPar)];                 
maintitle=myMainTitle;
savefile=sprintf('plots/E0_%s_qUScan_fixpar%s_%s.pdf',RunList,FixParLabel,mode);
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
hab=errorbar(myindex+2,ab_e0_stat,ab_e0err_stat,'o','Color',rgb('IndianRed'),'MarkerSize',7,'MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
hold on
hcw=errorbar(myindex-2,cw_e0_stat,cw_e0err_stat,'s','Color',rgb('DarkBlue'),'MarkerSize',7,'MarkerEdgeColor',rgb('DarkBlue'),'MarkerFaceColor',rgb('DarkBlue'),'LineWidth',2);
hold off
ylabel('Effective E_0 (eV)')
xlabel('Lower qU (V)'); 
leg=legend([hab,hcw],'Abdurashitov','CW_GLT','Location','NorthWest');                   
leg.Color = 'none'; leg.FontSize = 18; legend boxoff;
grid on
PrettyFigureFormat
set(gca,'FontSize',24);
grid on
xlim([-650 -50]);
publish_figure(gcf,'tmp.eps');