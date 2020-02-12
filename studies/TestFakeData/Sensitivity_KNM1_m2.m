% Quick Neutrino Mass Sensitivity Study for KNM1
% Compute and Display
CL=1;
RecomputeFlag = 'ON';
%Mode = 'FakeData';
Mode = 'Sim';
BkgFluct = 'Poisson';

% ASIMOV
if strcmp(RecomputeFlag,'ON')
%M=MultiRunAnalysis('DataType','Fake','FakeRunType','Fake1','RunList','KNM1_30d','exclDataStart',1,'fixPar','5 6','chi2','chi2Stat','Debug','ON');
M=MultiRunAnalysis('DataType','Fake','FakeRunType','Fake5y','RunList','KNM5y','exclDataStart',1,'fixPar','5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov');
M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;
end
%

% Fake Data
%M=MultiRunAnalysis('DataType','Fake','FakeRunType','Fake1','RunList','KNM1_30d','exclDataStart',1,'fixPar','5 6','chi2','chi2CMShape','Debug','ON');
%

chi2name = 'chi2CM';
exclDataStart = 1;
switch exclDataStart
    case 1
        belowE0 = 40;
end
switch Mode
    case 'Sim'
        savename = sprintf('./results/KNM1_NuMassSensitivity_%.0feV_%s.mat',belowE0,chi2name);
    case 'FakeData'
        savename = sprintf('./results/KNM1_NuMassSysBreakdownData_%.0feV_%s.mat',belowE0,chi2name);
end
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    StackCM = 'OFF';
    switch Mode
        case {'Sim'}
            M.RunData.qU     = M.ModelObj.qU;
            M.RunData.TBDIS  = M.ModelObj.TBDIS;
            M.RunData.TBDISE = M.ModelObj.TBDISE;
            StackCM = 'ON';
        case {'FakeData'}
            StackCM = 'ON';
    end
    M.fixPar = '5 6';
    M.exclDataStart = exclDataStart;
    Mass2      = zeros(8,1); % Single Sys: Stat, TC, TASR, FSD, RF
    Mass2Multi = zeros(8,1); % Multi: Stat, Stat+TC+TASR, +FSD, +RF
    %% Stat Sensitivity on Mass2:
    M.chi2 = 'chi2Stat';
    M.Fit;
    Mass2(1)      = M.FitResult.err(1);
    Mass2Multi(1) = M.FitResult.err(1);
    
    %Mass2(1)      = 0.3449;
    %Mass2Multi(1) = 0.3449;

    %% Systematics - does not work like that
%     M.FitCM_Obj.DTFlag = 'ON';
%     M.FitCM_Obj.TTFlag = 'ON';
%     M.FitCM_Obj.HTFlag = 'ON';
%     M.FitCM_Obj.WGTS_TASR_RelErr = 1e-3;
    
    %% With Single Systematic Effects
    M.chi2 = chi2name;
    
    fprintf(2,'Fit: Stat + TCoff_Rad + TCoff_OTHER + Stack \n');
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON'),'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2(2) = M.FitResult.err(1);
    Mass2Multi(2) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + TASR + Stack \n');
    M.ComputeCM('SysEffect',struct('TASR','ON'),'WGTS_TASR_RelErr',2e-4,'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2(3) =  M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + FSD + Stack \n');
    M.ComputeCM('SysEffect',struct('FSD','ON'),'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2(4) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + RF_BF + Stack \n');
    M.ComputeCM('SysEffect',struct('RF_BF','ON'),'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2(5) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + RF_EL + Stack \n');
    M.ComputeCM('SysEffect',struct('RF_EL','ON'),'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2(6) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + RF_RX + Stack \n');
    M.ComputeCM('SysEffect',struct('RF_RX','ON'),'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2(7) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + Bkg-Shape + Stack \n');
    M.ComputeCM('SysEffect',struct('TCoff_Rad','OFF','TCoff_OTHER','ON'),'Stack',StackCM,'BkgCM','ON','BkgFluct',BkgFluct);
    M.Fit;
    Mass2(8) = M.FitResult.err(1);
    
    %% With Stacked Systematics
    fprintf(2,'Fit: Stat + TCoff_Rad + TCoff_OTHER + TASR + Stack \n');
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON'),'WGTS_TASR_RelErr',2e-4,'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2Multi(3) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + TCoff_Rad + TCoff_OTHER + TASR + FSD + Stack \n');
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON','FSD','ON'),'WGTS_TASR_RelErr',2e-4,'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2Multi(4) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + TCoff_Rad + TCoff_OTHER + TASR + FSD + RF_BF + Stack \n');
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON','FSD','ON','RF_BF','ON'),'WGTS_TASR_RelErr',2e-4,'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2Multi(5) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + TCoff_Rad + TCoff_OTHER + TASR + FSD + RF_BF + RF_EL + Stack \n');
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON','FSD','ON','RF_BF','ON','RF_EL','ON'),'WGTS_TASR_RelErr',2e-4,'Stack',StackCM,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2Multi(6) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + TCoff_Rad + TCoff_OTHER + TASR + FSD + RF_RX + Stack \n');
    M.ComputeCM('Stack',StackCM,'WGTS_TASR_RelErr',2e-4,'BkgFluct',BkgFluct);
    M.Fit;
    Mass2Multi(7) = M.FitResult.err(1);
    
    fprintf(2,'Fit: Stat + TCoff_Rad + TCoff_OTHER + TASR + FSD + RF_RX + Stack + Bkg-Shape\n');
    M.ComputeCM('Stack',StackCM,'WGTS_TASR_RelErr',2e-4,'BkgCM','ON','BkgFluct',BkgFluct);
    M.Fit;
    Mass2Multi(8) = M.FitResult.err(1);
       
    %% for stacked plot
    Mass2Stacked = zeros(8,1);
    Mass2Single  = zeros(8,1);
    
    Mass2Stacked(1) = Mass2Multi(1);
    Mass2Single(1)  = Mass2(1);
    %Mass2Single(1)  = min(Mass2(:));
    
    for i=2:numel(Mass2Multi)
        Mass2Stacked(i) =  Mass2Multi(i)-Mass2Multi(i-1);
        Mass2Single(i)  =  Mass2(i)-Mass2(1);
    end
    StackEmpty = zeros(8,1);
    save(savename,'Mass2','Mass2Multi','Mass2Stacked','Mass2Single','StackEmpty');
end
%% Display
f55 = figure('Name','MultiBarPlot','Renderer','opengl');
set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
b = barh([10;20],1.00*[Mass2Stacked , StackEmpty]','stacked');
b(1).LineStyle = '--';
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';b(5).LineStyle = 'none';
b(1).FaceColor = rgb('White');
b(1).FaceColor = rgb('Amethyst');
b(3).FaceColor = rgb('LightCoral');
b(4).FaceColor = rgb('GoldenRod');
b(5).FaceColor = rgb('CadetBlue');
b(6).FaceColor = rgb('Navy');
b(7).FaceColor = rgb('PowderBlue');
b(8).FaceColor = rgb('DarkGreen');
leg_str = {' Statistical';'+ TC';'+ Tritium Acitivity Fluctuation';'+ Final State Distribution';...
    '+ Magnetic Fields';'+ Energy Loss Function'; '+ Column Density + Inel. Cross Section';'+ Background-Shape'}; %;'+ Response Function'};
b(1).BarWidth = b(1).BarWidth/2;
% Single Sys
SingleSysBarWidth = b(1).BarWidth/(6.4);
hold on;
pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
bTC  = barh([7.97-0.4;21], 1.00*[Mass2Single(1),Mass2Single(2); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bTC(2).FaceColor =rgb('Amethyst'); bTC(2).LineStyle = 'none'; bTC(2).FaceAlpha = 0.5;
bTC(1).FaceColor = 'w'; bTC(1).LineStyle = 'none'; %bFSD(1).FaceColor = 'w'; bRF(1).FaceColor = 'w';

bTASR  = barh([7.55-0.85;21], 1.00*[Mass2Single(1),Mass2Single(3); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bTASR(2).FaceColor =rgb('LightCoral'); bTASR(2).LineStyle = 'none'; bTASR(2).FaceAlpha = 0.5;
bTASR(1).FaceColor = 'w'; bTASR(1).LineStyle = 'none';

bFSD  = barh([6.62-0.85;21], 1.00*[Mass2Single(1),Mass2Single(4); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bFSD(2).FaceColor =rgb('GoldenRod'); bFSD(2).LineStyle = 'none'; bFSD(2).FaceAlpha = 0.5;
bFSD(1).FaceColor = 'w'; bFSD(1).LineStyle = 'none';

bRF_BF  = barh([5.63-0.85;21], 1.00*[Mass2Single(1),Mass2Single(5); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bRF_BF(2).FaceColor =rgb('CadetBlue'); bRF_BF(2).LineStyle = 'none'; bRF_BF(2).FaceAlpha = 0.5;
bRF_BF(1).FaceColor = 'w'; bRF_BF(1).LineStyle = 'none';

bRF_EL  = barh([4.57-0.85;21], 1.00*[Mass2Single(1),Mass2Single(6); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bRF_EL(2).FaceColor =rgb('PowderBlue'); bRF_EL(2).LineStyle = 'none'; bRF_EL(2).FaceAlpha = 0.9;
bRF_EL(1).FaceColor = 'w'; bRF_EL(1).LineStyle = 'none';

bRF_RX  = barh([3.44-0.85;21], 1.00*[Mass2Single(1),Mass2Single(7); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bRF_RX(2).FaceColor =rgb('DarkSlateGray'); bRF_RX(2).LineStyle = 'none'; bRF_RX(2).FaceAlpha = 0.5;
bRF_RX(1).FaceColor = 'w'; bRF_RX(1).LineStyle = 'none';

bBkgS  = barh([2.3-0.85;21], 1.00*[Mass2Single(1),Mass2Single(8); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bBkgS(2).FaceColor =rgb('DarkGreen'); bBkgS(2).LineStyle = 'none'; bBkgS(2).FaceAlpha = 0.5;
bBkgS(1).FaceColor = 'w'; bBkgS(1).LineStyle = 'none';

% Plot again for nicer statistical border --
b = barh([10;20],1.00*[Mass2Stacked , StackEmpty]','stacked');
b(1).LineStyle = 'none';b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';b(5).LineStyle = 'none'; b(6).LineStyle = 'none';b(7).LineStyle = 'none';b(8).LineStyle = 'none';
b(1).FaceColor = rgb('White');
b(2).FaceColor = rgb('Amethyst');
b(3).FaceColor = rgb('LightCoral');
b(4).FaceColor = rgb('GoldenRod');
b(5).FaceColor = rgb('CadetBlue');
b(6).FaceColor = rgb('PowderBlue');
b(7).FaceColor = rgb('DarkSlateGray');
b(8).FaceColor = rgb('DarkGreen');
b(1).BarWidth = b(1).BarWidth/2;

leg_str = {leg_str{:},'Single Contributions','Stat + Theoretical Corrections','Stat + Tritium Activity Fluctuation', 'Stat + Final State Distribution',...
    'Stat + Magnetic Fields','Stat + Energy Loss Function','Stat + Column Density + Inel. Cross Section','Stat + Background Shape'};
leg = legend([b(1),b(2),b(3),b(4),b(5),b(6),b(7),b(8),pnone,bTC(2),bTASR(2),bFSD(2),bRF_BF(2),bRF_EL(2),bRF_RX(2),bBkgS(2)],leg_str{:});
leg.Location = 'north';
leg.NumColumns = 2; legend boxoff
PrettyFigureFormat;
ylim([0. 18]);
%xlim(1.00*[0.337,0.4]);
%xlim(1.00*[0.45,.55]);
xlim(1.00*[0.06,.115]);
%xlim(1.00*[0.1,.2]);

grid on;
set(gca,'FontSize',18);
switch Mode
    case 'Sim'
%        xlabel(sprintf('m^2 - 68.3%% C.L. Upper Limit assuming m^2 = 0 eV^2/c^2'));
        xlabel(sprintf('m^2 (eV/c^2)^2 - 68.3%% C.L. Uncertainty'));
    case 'FakeData'
        xlabel(sprintf('m^2 (eV/c^2)^2 - 68.3%% C.L. Uncertainty'));
end
yticks = [];yticklabels('');
%title('KNM1 30d - 124 Pixels Stacked - 360 runs Stacked -  neutrino mass squared sensitivity Breakdown'); 
title('KATRIN 900 days - 124 Pixels Stacked - 90 runs Stacked -  neutrino mass squared sensitivity Breakdown'); 

% Text - 1 sigma error on mass squared +-stat +- sys
mstat = Mass2Multi(1);
msys  = sqrt(Mass2Multi(8)^2-Mass2Multi(1)^2);
strms = sprintf('m^2 = 0 \\pm %0.2f (stat) \\pm %0.2f (sys) eV^2 \n',mstat,msys);
%text(0.48,12.,strms,'FontSize',20,'FontWeight','bold','FontName','Arial');
%text(0.71,12.,strms,'FontSize',20,'FontWeight','bold','FontName','Arial');
text(0.08,12.,strms,'FontSize',20,'FontWeight','bold','FontName','Arial');

%% Save
unix('mkdir ./plots/png')
plot_name = sprintf('./plots/png/KNM1_Mass2SensitivityBreakdown_%s.png',Mode);
print(f55,plot_name,'-dpng');

%%




