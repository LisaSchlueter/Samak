% Compute the probability to reject the no electronics excited states
% hypothesis based on a frequentist montecarlo simulation
% for a given list of First Tritium data runs
%   
% Call RunWiseFSDeeonsetRejection.m
% 
% T. Lasserre, Pablo Morales, Lisa Schlueter
% Last Update: June 24 2018
%

exclDataStart = 7;
chi2          = 'stat';
index =  [1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];

[ n, datafiles , RunList ] = GetRunList( '../../tritium-data/hdf5/','*.h5',1,'string');
RunList=RunList(RunList>40666);
RunList=RunList(RunList<40693);

r=[];PnoEE=zeros(1,numel(RunList));
pvaluenoEE=zeros(1,numel(RunList));
pvalueEE=zeros(1,numel(RunList));
for i=1:1:numel(RunList)
    fprintf(2,'processing run %g ...\n',RunList(i))
    [PnoEE(i) , fit_p_data_ee , fit_p_data_noee]= RunWiseFSDeeonsetRejection('run',RunList(i),'exclDataStart',exclDataStart,'runSim','OFF');
    pvaluenoEE(i) = chi2pvalue(fit_p_data_noee(13),fit_p_data_noee(14));
    pvalueEE(i) = chi2pvalue(fit_p_data_ee(13),fit_p_data_ee(14));
    fprintf(2,'FSD No EE hypothesis rejected at  %g %% CL ...\n',PnoEE(i))
    r    = [r num2str(RunList(i))]; r    = string(r);
end

%% Display (rejection probability) 
maintitle=sprintf('KATRIN  First Tritium - All pixels but 2 ext Rings - FSD EE Onset');
savefile=sprintf('plots/KATRIN_FT_AllPixels_Samak_%s_%geVbelowE_0_FSDonset_Pnoee.pdf',chi2,index(exclDataStart));
x   = linspace(1,numel(RunList),numel(RunList));
gcf = figure('Name','VFT FSD Onset','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
errorbar(1:1:numel(RunList),PnoEE,PnoEE.*0,'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
xticks(x)
xticklabels(r)
ylabel('Rejection Probability (%)')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
ylim([0. 110])
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
publish_figurePDF(gcf,savefile);

%% Display (pvalues)   
aintitle=sprintf('KATRIN  First Tritium - All pixels but 2 ext Rings - FSD EE Onset');
savefile=sprintf('figures/KATRIN_FT_AllPixels_Samak_%s_%geVbelowE_0_FSDonset_Pnoee2.pdf',chi2,index(exclDataStart));
x   = linspace(1,numel(RunList),numel(RunList));
gcf = figure('Name','VFT FSD Onset','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
hee=errorbar(1:1:numel(RunList),pvalueEE*100,pvalueEE.*0,'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
hold on
hnoee=errorbar(1:1:numel(RunList),pvaluenoEE*100,pvaluenoEE.*0,'ks','MarkerSize',12,'MarkerFaceColor',rgb('Khaki'),'LineWidth',1);
hold off
xticks(x)
xticklabels(r)
ylabel('p-value (%)')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
ylim([0. 110])
grid on
legend([hnoee hee],'p-value fit without Electronic Excited States','p-value fit with Electronic Excited States');
set(gca,'yscale','log');
set(gca,'FontSize',16);
PrettyFigureFormat
publish_figurePDF(gcf,savefile);
