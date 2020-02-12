%  Fit the probability to decay to the HeD ground state 
%  Call FT_FitFSD for a set of First Tritium runs
%  (obtained via GetRunList, see code)
%  - E0 , B , N free, m=0, Pgs free
%
% Input:
%  - chi2            stat, 
%  - exclDataStart:  7=400eV, 9=200eV, ... 
% Output:
%  - matfile including:
%    fit_p_data_ee: fit parameters
%    Pgs: fit probabilities Pgs
%    PgsErr: fit probabilities Pgs errors 
% 
% T. Lasserre, Pablo Morales, Lisa Schlueter
% Last Update: June 24 2018
%

exclDataStart = 7;
chi2          = 'stat';
index =  [1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];
fit   = 'ON';

[ n, datafiles , RunList ] = GetRunList( '../../tritium-data/hdf5/','*.h5',1,'string');
RunList=RunList(RunList>40666);
RunList=RunList(RunList<40693);

Sref = ref_runsummaries(RunList(1),'ex2','ISCS','Theory','recomputeRF','OFF');
switch fit
    case 'ON'
r=[];Pgs=zeros(1,numel(RunList)); PgsErr=zeros(1,numel(RunList));
for i=1:1:numel(RunList)
    fprintf(2,'processing run %g ...\n',RunList(i));
    [fit_p_data_ee(i,:) , ~]  =  FT_FitFSD('run',RunList(i),'FSD','eeFit','Mode','Data','nfit',1,'exclDataStart',exclDataStart,'chi2','chi2Stat');
    Pgs(i) = Sref.DTNormGS_i+fit_p_data_ee(i,5); PgsErr(i) = fit_p_data_ee(i,11);
    fprintf(2,'FSD Onset : Probability to decay to GS = %g +- %g ...\n',Pgs(i),PgsErr(i))
    r    = [r num2str(RunList(i))]; r    = string(r);
end
end

%% Display
maintitle=sprintf('KATRIN  First Tritium - All pixels - FSD <P_{gs,DT}>=%.2f \\pm %.2f (std)',mean(Pgs),std(Pgs));
savefile=sprintf('plots/KATRIN_FT_AllPixels_Samak_%s_%geVbelowE_0_FSDonset_Pgs_1.pdf',chi2,index(exclDataStart));
x   = linspace(1,numel(RunList),numel(RunList));
gcf = figure('Name','VFT FSD Onset','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
errorbar(1:1:numel(RunList),Sref.DTNormGS_i+fit_p_data_ee(:,5),fit_p_data_ee(:,11),'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
hold on
line([min(x)-0.5,max(x)+0.5],[Sref.DTNormGS_i Sref.DTNormGS_i],'Color','Red','LineWidth',1,'LineStyle','--')
hold off
xticks(x)
xticklabels(r)
ylabel('Ground State Decay Probability (%)')
xlabel('run');xtickangle(45);
xlim([0.5 numel(RunList)+0.5])
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
publish_figurePDF(gcf,savefile);

%% Analysis - Pgs verus Chi2
maintitle=sprintf('KATRIN  First Tritium - All pixels - FSD <P_{gs,DT}>=%.2f \\pm %.2f (std)',mean(Pgs),std(Pgs));
savefile=sprintf('plots/KATRIN_FT_AllPixels_Samak_%s_%geVbelowE_0_FSDonset_Pgs_2.pdf',chi2,index(exclDataStart));
gcf = figure('Name','VFT FSD Onset','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
errorbar(Sref.DTNormGS_i+fit_p_data_ee(:,5),fit_p_data_ee(:,13),zeros(numel(RunList),1),'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
xlabel('Ground State Decay Probability (%)')
ylabel('\chi2');xtickangle(45);
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
publish_figurePDF(gcf,savefile);

%% Analysis - Pgs verus E0
maintitle=sprintf('KATRIN  First Tritium - All pixels - FSD <P_{gs,DT}>=%.2f \\pm %.2f (std)',mean(Pgs),std(Pgs));
savefile=sprintf('figures/KATRIN_FT_AllPixels_Samak_%s_%geVbelowE_0_FSDonset_Pgs_3.pdf',chi2,index(exclDataStart));
gcf = figure('Name','VFT FSD Onset','NumberTitle','off','rend','painters','pos',[10 10 1200 600]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
errorbar(18573.7+fit_p_data_ee(:,2),Sref.DTNormGS_i+fit_p_data_ee(:,5),fit_p_data_ee(:,5).*0,'ks','MarkerSize',12,'MarkerFaceColor',rgb('CadetBlue'),'LineWidth',1);
ylabel('Ground State Decay Probability (%)')
xlabel('(E_0)_{eff}-18575');xtickangle(45);
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
publish_figurePDF(gcf,savefile);

%% Savez data file
switch fit
    case 'ON'
savefile=sprintf('figures/KATRIN_FT_AllPixels_Samak_%s_%geVbelowE_0_FSDonset_Pgs_NoE0Pull.mat',chi2,index(exclDataStart));
save(savefile,'fit_p_data_ee','Pgs','PgsErr');
end
