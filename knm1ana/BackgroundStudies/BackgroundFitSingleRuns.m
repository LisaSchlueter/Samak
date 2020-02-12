%
% Macro to Study Quality of Run-Wise Fit in Background Region
% T. Lasserre
% Last Modified: 4/6/19
%

counter = 0;

if exist('KNM1_runWiseBackgroundFits.mat','file')
    sprintf('FitResults Stat Only already computed - Loading from File \n')
    load('KNM1_runWiseBackgroundFits.mat');
else
    
    for i=Real.RunList'
        
        counter  = counter + 1;
        
        A=RunAnalysis('RunNr',i,'DataType','Real',...
            'fixPar','1 5 6 7 8 9 10 11',...
            'i_qUOffset',0,'exclDataStart',14,...
            'NonPoissonScaleFactor',1,...
            'fitter','matlab');
        A.Fit ;
        %A.PlotFit('saveplot','OFF') ;
        
        A=RunAnalysis('RunNr',i,'DataType','Real',...
            'fixPar','1 2 3 4 5 6 7 8 9 10 11',...
            'i_mnu',0,...
            'i_Q',A.FitResult.par(2),...
            'i_B',A.FitResult.par(3),...
            'i_N',A.FitResult.par(4),...
            'i_qUOffset',0,...
            'exclDataStart',36,...
            'NonPoissonScaleFactor',1,...
            'fitter','matlab');
        A.Fit ;
        %A.PlotFit('saveplot','OFF') ;
       
        % Save Results
        Results(counter) = A.FitResult;
    end
    save('KNM1_runWiseBackgroundFits.mat','Real','Results');
end

index     = A.ModelObj.Q_i-A.ModelObj.qU;
Results_T = struct2table(Results);
pv        = 1-chi2cdf(Results_T.chi2min,Results_T.dof);

%% Chi2 distribution
myMainTitle=[sprintf('KATRIN - Run-wise Background Fits - %d Runs (%s - %.1f',...
    numel(Real.RunList),A.chi2,abs(double(index(A.exclDataStart)))),'eV below endpoint)'];
maintitle=myMainTitle;
savefile=sprintf('plots/MRA_RunWiseBackgroundFitResults%d_%s_%.0feVbelowE0-1.png',...
    numel(Real.RunList),A.chi2,abs(index(A.exclDataStart)));
fig1 = figure('Name','MRA Background Run-wise Fits','NumberTitle','off','rend','painters','pos',[10 10 1000 600]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';

e=histogram(chi2rnd(5,1e6,1),'normalization', 'pdf','FaceAlpha',0.5,'FaceColor',rgb('SteelBlue'));
hold on
d=histogram(Results_T.chi2min,'normalization', 'pdf','FaceAlpha',0.5,'FaceColor',rgb('IndianRed'));
hold off
xlabel('\chi ^2');
ylabel('pdf');
legend([e d],'Expected from \chi^2 statistics','KNM1 background data - 274 runs','Location','northeast');
PrettyFigureFormat
set(gca,'yscale','log');
set(gca,'FontSize',18);
xlim([-0.5 30]);
export_fig(gcf,savefile,'-m3');

%% Chi2 distribution - Versus Run
% Labels
x   = linspace(1,numel(Real.RunList),numel(Real.RunList));
index = Real.ModelObj.Q_i-Real.ModelObj.qU;
r='';
for i=1:numel(Real.RunList)
    r    = [r num2str(Real.RunList(i))];
    r    = string(r);
end


myMainTitle=[sprintf('KATRIN - Run-wise Background Fits - %d Runs (%s - %.1f',...
    numel(Real.RunList),A.chi2,abs(double(index(A.exclDataStart)))),'eV below endpoint)'];
maintitle=myMainTitle;
savefile=sprintf('plots/MRA_RunWiseBackgroundFitResults%d_%s_%.0feVbelowE0-2.png',...
    numel(Real.RunList),A.chi2,abs(index(A.exclDataStart)));
fig1 = figure('Name','MRA Background Run-wise Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 2000]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');
a.FontSize=20;a.FontWeight='bold';
                   
stairs(x,pv,'s','MarkerSize',10,'MarkerEdgeColor',rgb('IndianRed'),'MarkerFaceColor',rgb('IndianRed'),'LineWidth',1)
set(gca,'yscale','log');
xticks(x(1:5:end)); xticklabels(r(1:5:end));
xtickangle(45);
MyMarkerSize=8;
PrettyFigureFormat
export_fig(gcf,savefile,'-m3');

%% Runs with largest chi2
[a b] = sort(Results_T.chi2min);
TagRuns = Real.RunList(b);
TagChi2 = a;
for i=1:1:274
fprintf('Run: %s - chi2 = %0.1f (6 data points) \n',num2str(TagRuns(i)),TagChi2(i));
end