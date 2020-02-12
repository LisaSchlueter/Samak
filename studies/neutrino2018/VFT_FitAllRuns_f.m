function  VFT_FitAllRuns_f(varargin)

close all;

addpath(genpath('../../../Samak2.0'));

    % Parser
    p = inputParser;
    p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('exclDataStart',7,@(x)isfloat(x) && x>=0);
    p.addParameter('chi2','CM',@(x)ismember(x,{'CM','Stat'}));
    p.parse(varargin{:});
    display           =    p.Results.display;
    exclDataStart     =    p.Results.exclDataStart; % Start Fit at bin i (i=1 -> qU(min)=-1.7keV, i=9 ->qU(min)=200eV, i=7 ->qU(min)=400eV)
    chi2              =    p.Results.chi2;
    
%
index =  [1602  1402  1204  1002   804   602   402   304   202   177   154   129   104    92    82    72    62    52    45    32    22    15     2    -6   -18   -28];

%
[par40257 , err40257 , chi2min40257 , dof40257] = VFT_FitRunNNu18_f('display','OFF','RunNr',40257,'exclDataStart',exclDataStart,'chi2',chi2);
[par40258 , err40258 , chi2min40258 , dof40258] = VFT_FitRunNNu18_f('display','OFF','RunNr',40258,'exclDataStart',exclDataStart,'chi2',chi2);
[par40259 , err40259 , chi2min40259 , dof40259] = VFT_FitRunNNu18_f('display','OFF','RunNr',40259,'exclDataStart',exclDataStart,'chi2',chi2);
[par40260 , err40260 , chi2min40260 , dof40260] = VFT_FitRunNNu18_f('display','OFF','RunNr',40260,'exclDataStart',exclDataStart,'chi2',chi2);
[par40263 , err40263 , chi2min40263 , dof40263] = VFT_FitRunNNu18_f('display','OFF','RunNr',40263,'exclDataStart',exclDataStart,'chi2',chi2);
[par40264 , err40264 , chi2min40264 , dof40264] = VFT_FitRunNNu18_f('display','OFF','RunNr',40264,'exclDataStart',exclDataStart,'chi2',chi2);
[par40265 , err40265 , chi2min40265 , dof40265] = VFT_FitRunNNu18_f('display','OFF','RunNr',40265,'exclDataStart',exclDataStart,'chi2',chi2);
[par40266 , err40266 , chi2min40266 , dof40266] = VFT_FitRunNNu18_f('display','OFF','RunNr',40266,'exclDataStart',exclDataStart,'chi2',chi2);

%
% [par40257cm , err40257cm , chi2min40257cm , dof40257cm] = VFT_FitRunNNu18_f('display','OFF','RunNr',40257,'exclDataStart',7,'chi2','CM');
% [par40258cm , err40258cm , chi2min40258cm , dof40258cm] = VFT_FitRunNNu18_f('display','OFF','RunNr',40258,'exclDataStart',7,'chi2','CM');
% [par40259cm , err40259cm , chi2min40259cm , dof40259cm] = VFT_FitRunNNu18_f('display','OFF','RunNr',40259,'exclDataStart',7,'chi2','CM');
% [par40260cm , err40260cm , chi2min40260cm , dof40260cm] = VFT_FitRunNNu18_f('display','OFF','RunNr',40260,'exclDataStart',7,'chi2','CM');
% [par40263cm , err40263cm , chi2min40263cm , dof40263cm] = VFT_FitRunNNu18_f('display','OFF','RunNr',40263,'exclDataStart',7,'chi2','CM');
% [par40264cm , err40264cm , chi2min40264cm , dof40264cm] = VFT_FitRunNNu18_f('display','OFF','RunNr',40264,'exclDataStart',7,'chi2','CM');
% [par40265cm , err40265cm , chi2min40265cm , dof40265cm] = VFT_FitRunNNu18_f('display','OFF','RunNr',40265,'exclDataStart',7,'chi2','CM');
% [par40266cm , err40266cm , chi2min40266cm , dof40266cm] = VFT_FitRunNNu18_f('display','OFF','RunNr',40266,'exclDataStart',7,'chi2','CM');

x   = [1 2 3 4 5 6 7 8 ];
l   = ['40257' '40258' '40259' '40260' '40263' '40264' '40265' '40266'];

%% Stat
% Endpoint
E0    = [par40257(2) par40258(2) par40259(2) par40260(2) par40263(2) par40264(2) par40265(2) par40266(2)];
E0Err = [err40257(2) err40258(2) err40259(2) err40260(2) err40263(2) err40264(2) err40265(2) err40266(2)];
% Background
B    = [par40257(3) par40258(3) par40259(3) par40260(3) par40263(3) par40264(3) par40265(3) par40266(3)];
BErr = [err40257(3) err40258(3) err40259(3) err40260(3) err40263(3) err40264(3) err40265(3) err40266(3)];
% Normalizationn
N    = [par40257(4) par40258(4) par40259(4) par40260(4) par40263(4) par40264(4) par40265(4) par40266(4)];
NErr = [err40257(4) err40258(4) err40259(4) err40260(4) err40263(4) err40264(4) err40265(4) err40266(4)];
% Chi2
Chi2 = [chi2pvalue(chi2min40257,dof40257) chi2pvalue(chi2min40258,dof40258) chi2pvalue(chi2min40259,dof40259) chi2pvalue(chi2min40260,dof40260) chi2pvalue(chi2min40263,dof40263) chi2pvalue(chi2min40264,dof40264) chi2pvalue(chi2min40265,dof40265) chi2pvalue(chi2min40266,dof40266)];


%% Plot
switch display
    case 'ON'
maintitle=sprintf('KATRIN Very First Tritium - All pixels - Samak Fit: %s , %g eV below E_0',chi2,index(exclDataStart));
savefile=sprintf('figures/KATRIN_VFT_AllPixels_Samak_%s_%geVbelowE_0.eps',chi2,index(exclDataStart));
x   = [1 2 3 4 5 6 7 8 ];
r   = {'40257' '40258' '40259' '40260' '40263' '40264' '40265' '40266'};
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';
subplot(2,2,1)
bar(x,E0-18575,'facecolor',[1 0.549 0]); 
hold on
errorb(x,E0-18575,E0Err,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
hold off
xticks(x)
xticklabels(r)
ylabel('E_0-18575 (eV)')
xlabel('run')
xlim([0.5 8.5]) 
sp1 = sprintf('<E0>=%.2f eV \\pm %.2f eV (std)',mean(E0),std(E0-18575));title(sp1) 
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
subplot(2,2,2)
bar(x,B,'facecolor',[1 0.549 0]); 
hold on
errorb(x,B,BErr,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
hold off
xticks(x)
xticklabels(r)
ylabel('Background (cps)')
xlabel('run')
xlim([0.5 8.5]) 
sp2 = sprintf('<B>=%.2f cps \\pm %.2f cps (std)',mean(B),std(B));title(sp2) 
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
subplot(2,2,3)
bar(x,N,'facecolor',[1 0.549 0]); 
hold on
errorb(x,N,NErr,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
hold off
xticks(x)
xticklabels(r)
ylabel('1+Normalization')
xlabel('run')
xlim([0.5 8.5]) 
sp3 = sprintf('<1+N>=%.2f \\pm %.2f  (std)',mean(N),std(N));title(sp3) 
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
subplot(2,2,4)
bar(x,Chi2,'facecolor',[1 0.549 0])
hold on
line([min(x)-0.5,max(x)+0.5],[0.05 0.05],'Color','Red','LineWidth',2,'LineStyle','--')
hold off
xticks(x)
xticklabels(r)
ylabel('p-value')
xlabel('run')
xlim([0.5 8.5]) 
sp3 = sprintf('<p-value>=%.2f \\pm %.2f (std)',mean(Chi2),std(Chi2));title(sp3) 
set(gca,'yscale','log');
grid on
set(gca,'FontSize',16);
PrettyFigureFormat
publish_figure(gcf,savefile);
end

