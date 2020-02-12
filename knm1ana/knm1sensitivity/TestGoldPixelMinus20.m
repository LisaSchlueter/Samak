%
% Statitical Test
% KNM1 Analysis
% Remove randomly 20 pixels from Golden Run List
% Check probability for a Given Neutrino Mass Shift
%
% Thierry Lasserre
% Last Modified: 26/07/2019
%

%% settings
RunList               = 'KNM1';
exclDataStart         = 14;
nTrials               = 250 ;
RecomputeFlag         = 'OFF';
SysEffects            = struct('TASR','ON','FSD','ON','RF_RX','ON','RF_EL','ON','RF_BF','ON','BkgShape','ON','TCoff_RAD','ON','TCoff_OTHER','ON','Stack','ON');
BkgCM                 = 'ON';
%%
for i=1:1:nTrials
    
%% Init Model Object and covariance matrix object
Real = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2CMShape','DataType','Real',...
    'exclDataStart',exclDataStart,...
    'fixPar','5 6 7 8 9 10 11',...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','migrad');
Real.ComputeCM('SysEffects',SysEffects,'BkgCM',BkgCM);
Real.Fit('CATS','OFF'); % Real.PlotFit('saveplot','png');
Rsysshape(i) = Real.FitResult;

end
%save('TestGoldPixelMinus20.mat','Rsysshape');

%% Extracting and plotting the data
mass2=[];
mass2E=[];
for i=1:1:nTrials
tmp1 = Rsysshape(i).par;
tmp1 = tmp1(1);
mass2=[mass2 tmp1];

tmp2 = Rsysshape(i).err;
tmp2 = tmp2(1);
mass2E=[mass2E tmp2];

end

%% Mass Squared Shift
myMainTitle = sprintf('KATRIN KNM1 Fit including 97 random pixels of the Golden Run List');
maintitle   = myMainTitle;
savefile    = sprintf('TestGoldPixelMinus20_1.png');
fig1      = figure('Name','KATRIN KNM1 Fit including 97 random pixels of the Golden Run List','NumberTitle','off','rend','painters','pos',[10 10 1300 800]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

t=histogram(mass2,250,'Normalization','cdf','LineWidth',2)
hold on

r20ext = line([-0.257 -0.257],[0 1],'Color','Red','LineWidth',3);
[i ii] = find(t.BinEdges>=-0.257); l = t.Values(min(ii));
r20ext = line([t.BinEdges(1) t.BinEdges(end)],[l l],'LineStyle','--','Color','Red','LineWidth',3);

r = line([-0.95 -0.95],[0 1],'Color','Blue','LineWidth',3);
hold off
grid on
xlabel('m^2 Best Fit (eV^2)');
ylabel('Cumulative Probability');
%title('KNM1 Fit including 97 random pixels of the Golden Run List ');
legend([t,r,r20ext],'Trials','KNM1 Golden Run List','KNM1 with 9 inner rings');
PrettyFigureFormat
ylim([0 1]);
xlim([t.BinEdges(1) t.BinEdges(end)]);


%% Mass Squared Error
myMainTitle = sprintf('KATRIN KNM1 Fit including 97 random pixels of the Golden Run List');
maintitle   = myMainTitle;
savefile    = sprintf('TestGoldPixelMinus20_1.png');
fig1      = figure('Name','KATRIN KNM1 Fit including 97 random pixels of the Golden Run List','NumberTitle','off','rend','painters','pos',[10 10 1300 800]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

t=histogram(mass2E,250,'Normalization','cdf','LineWidth',2)
hold on

r20ext = line([0.94 0.94],[0 1],'Color','Red','LineWidth',3);
[i ii] = find(t.BinEdges>=0.94); l = t.Values(min(ii));
r20ext = line([t.BinEdges(1) t.BinEdges(end)],[l l],'LineStyle','--','Color','Red','LineWidth',3);

r = line([1.02 1.02],[0 1],'Color','Blue','LineWidth',3);
hold off
grid on
xlabel('m^2 Error from Fit (eV^2)');
ylabel('Cumulative Probability');
%title('KNM1 Fit including 97 random pixels of the Golden Run List ');
legend([t,r,r20ext],'Trials','KNM1 Golden Run List','KNM1 with 9 inner rings');
PrettyFigureFormat
ylim([0 1]);
xlim([t.BinEdges(1) t.BinEdges(end)]);


%% Mass Squared Error / Mass2 2D
myMainTitle = sprintf('KATRIN KNM1 Fit including 97 random pixels of the Golden Run List');
maintitle   = myMainTitle;
savefile    = sprintf('TestGoldPixelMinus20_1.png');
fig1      = figure('Name','KATRIN KNM1 Fit including 97 random pixels of the Golden Run List','NumberTitle','off','rend','painters','pos',[10 10 1300 800]);
a=annotation('textbox', [0 0.9 1 0.1], 'String', maintitle,'EdgeColor', 'none','HorizontalAlignment', 'center');
a.FontSize=24;a.FontWeight='bold';

t=scatter(mass2,mass2E,...
    'MarkerEdgeColor',rgb('SteelBlue'),...
    'MarkerFaceColor',rgb('CadetBlue'),...
    'LineWidth',2)

r20ext = line([-2.5 2.5],[0.94 0.94],'Color','Red','LineWidth',3,'LineStyle','--');
r20ext = line([-0.257 -0.257],[-1 3],'Color','Red','LineWidth',3,'LineStyle','--');
          grid on
ylabel('m^2 Error from Fit (eV^2)');
xlabel('m^2 Best Fit (eV^2)');
PrettyFigureFormat
ylim([0.7 1.4]);
xlim([-2.4 0.1]);
