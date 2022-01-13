% Diagnosis plots for energy loss function covariance matrix

%% settings
RunList = 'KNM1';%'FTpaper';
ELossFlag ='KatrinT2'; %'Abdurashitov';;%'KatrinD2';%'Abdurashitov';%''Aseev';%Abdurashitov';%KatrinD2';%'Abdurashitov';%'KatrinD2';
nTrials = 1000;
RecomputeFlag = 'OFF';
Compare = 'OFF'; %compare with fitrium errorband (abdurashitov)
nScat = 1; % 1 or 2
% Init Model Object and covariance matrix object
if ~exist('A','var')
    A = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','NonPoissonScaleFactor',1);
end

A.ModelObj.ELossFlag = ELossFlag; A.ModelObj.InitializeELossFunction;
SysEffects = struct('RF_EL','ON');
CM = CovarianceMatrix('StudyObject',A.ModelObj, 'nTrials',nTrials,'SysEffect',SysEffects,'RecomputeFlag',RecomputeFlag);

if strcmp(Compare,'ON') && strcmp(ELossFlag,'Abdurashitov')
    d = importdata('Fitrium_ELossAbdurashitovBand.dat');
    E_f         = d.data(:,1); %energy fitrium
    Eloss_f     = d.data(:,2); %eloss fitrium
    ElossErr_f  = d.data(:,3); %eloss std fitrium
    maxE = 50;%round(max(E_f),0);   
    xmin = 8; xmax = 20;
    E = E_f;
else
    maxE = 200;
    xmin = 5; xmax = 25;
    minE=-maxE; NbinE = (maxE-minE)/0.04;
    E = linspace(minE,maxE,NbinE);
end
Estep = E(2) - E(1);
[~,fscatnE] = CM.ComputeELossLookupTable('E',E);
%% plots first two scatterings
f33 = figure('Renderer','opengl');
set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);

% just take short range for plot

plotE =E(E>xmin & E<xmax);
Eloss1 = squeeze(fscatnE(1,E>xmin & E<xmax,:));
Eloss2 = squeeze(fscatnE(2,E>xmin & E<xmax,:));

if ismember(ELossFlag,{'KatrinD2','KatrinT2'})
    std1 = std(Eloss1')*10;
    std2 = std(Eloss2')*10;
    legstr = '1\sigma error band x 10';
else
    std1 = std(Eloss1');
    std2 = std(Eloss2');
    legstr = '1\sigma error band';
  
end

% plot
[l1, a1] = boundedline(plotE,mean(Eloss1,2),std1); l1.Color = rgb('FireBrick'); a1.FaceColor = rgb('FireBrick'); a1.FaceAlpha = 0.2;
l1.LineWidth = 2;
hold on;
 ymax = max(mean(Eloss1,2))*1.1;
if strcmp(Compare,'ON')
    [l2, a2] = boundedline(E_f,Eloss_f,ElossErr_f); l2.Color = rgb('SteelBlue'); a2.FaceColor = rgb('SteelBlue'); a2.FaceAlpha = 0.2; l2.LineWidth = 2;
    l2.LineStyle = '--';
    legstr12 = 'Samak';
    legstr22 = 'Fitrium';
    mylocation = 'northeast';
    
    %[~,fscatnE]= A.ModelObj.ComputeELossFunction('E',E_f);
    
    %     maxE = 9288;
    %     minE=-maxE; NbinE = (maxE-minE)/0.1;
    %     Ep = linspace(minE,maxE,NbinE);
    %plot(E_f,fscatnE(1,:),'LineWidth',3,'Color','k');
else
    legstr12 = '1st  scattering';
    if nScat==2
        [l2, a2] = boundedline(plotE,mean(Eloss2,2),std2); l2.Color = rgb('SteelBlue'); a2.FaceColor = rgb('SteelBlue'); a2.FaceAlpha = 0.2; l2.LineWidth = 2;
        legstr22 = '2nd  scattering';
        ymax = max(mean(Eloss2,2))*1.1;
    end
end

if nScat==2 || strcmp(Compare,'ON')
    mylocation = 'northwest';
    leg = legend([l1,l2,a1,a2],legstr12, legstr22,legstr,legstr,...
        'Location', mylocation); legend boxoff
    leg.NumColumns=2;
else
   mylocation = 'northwest';
    leg = legend([l1,a1],legstr12,legstr,...
        'Location', mylocation); legend boxoff  
end
PrettyFigureFormat;
xlim([xmin xmax]);
ylim([0 ymax])
xlabel(['energy loss ',char(949),' (eV)']);
ylabel('probability density (eV^{-1})')
title(sprintf('energy loss function %s with uncertainty band',A.ModelObj.ELossFlag));


% save plot
savepath = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
switch Compare
    case 'OFF'
savename = [savepath,sprintf('knm1_eloss_%s.png',A.ModelObj.ELossFlag)];
    case 'ON'
savename = [savepath,sprintf('knm1_eloss_%s_Comparison_%.0f.png',A.ModelObj.ELossFlag,nScat)];
end
if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end
print(gcf,savename,'-dpng','-r400');

