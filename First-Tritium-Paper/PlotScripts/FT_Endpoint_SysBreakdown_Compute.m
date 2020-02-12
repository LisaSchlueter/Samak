% Quick Neutrino E0 Sensitivity Study for First Tritium
% Compute and Display
Mode = 'Sim';
RecomputeFlag = 'OFF';
%Mode = 'Data';
RunList = 'FTpaper';
exclDataStart = 13;
if ~exist('M','var')
M = MultiRunAnalysis('RunList',RunList,'exclDataStart',exclDataStart,'DataType','Real',...
    'FSDFlag','SAENZ','ELossFlag','Abdurashitov','SysBudget',0,'DataEffCor','RunSummary',...
    'NonPoissonScaleFactor',1);
end
%% Fit);
M.ModelObj.ComputeTBDDS; M.ModelObj.ComputeTBDIS;

chi2name = 'chi2CMShape';
switch exclDataStart
    case 1
        belowE0 = 1600;
    case 7 
        belowE0 = 400;
    case 9
        belowE0 = 200;
    case 12
        belowE0 = 125;
    case 13
        belowE0 = 100;
end
%mydir = [getenv('SamakPath'),'studies/local_LisaThesisPlots/results/'];
mydir = [getenv('SamakPath'),'First-Tritium-Paper/PlotScripts/results/'];
if ~exist(mydir,'dir')
    system(['mkdir ',mydir]);
end

switch Mode
    case 'Sim'
        savename = [mydir,sprintf('FT_EndpointSensitivity_%.0feV_%s_%s.mat',belowE0,chi2name,RunList)];
    case 'Data'
        savename = [mydir,sprintf('FT_EndpointSysBreakdownData_%.0feV_%s_%s.mat',belowE0,chi2name,RunList)];
end
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename)
else
    StackCM = 'OFF';
    switch Mode
        case 'Sim'
            M.RunData.qU = M.ModelObj.qU;
            M.RunData.TBDIS = M.ModelObj.TBDIS;
            M.RunData.TBDISE = M.ModelObj.TBDISE;
            StackCM = 'OFF';
    end
    M.fixPar = '1 5 6 7 8 9 10 11';
    M.exclDataStart = exclDataStart;
    E0      = zeros(9,1); % Single Sys: Stat, TC, TASR, FSD, BF,EL,RX,Stack, Bkg
    E0Multi = zeros(9,1); % Multi: Stat, Stat+TC, +TASR, +FSD, +BF, +EL, +RX,+Stack, +bkg
    %% Stat Sensitivity on E0:
    M.Fit;
    E0(1) = M.FitResult.err(2);
    E0Multi(1) = M.FitResult.err(2);
    
    %% With Single Systematic Effects
    M.chi2 = chi2name;
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON'),'BkgCM','OFF');
    M.Fit;
    E0(2) = M.FitResult.err(2);
    E0Multi(2) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('TASR','ON'),'BkgCM','OFF');
    M.Fit;
    E0(3) =  M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('FSD','ON'),'BkgCM','OFF');
    M.Fit;
    E0(4) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('RF_BF','ON'),'BkgCM','OFF');
    M.Fit;
    E0(5) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('RF_EL','ON'),'BkgCM','OFF');
    M.Fit;
    E0(6) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('RF_RX','ON'),'BkgCM','OFF');
    M.Fit;
    E0(7) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('FPDeff','ON'),'BkgCM','OFF');
    M.Fit;
    E0(8) = M.FitResult.err(2); 
    
    M.ComputeCM('SysEffect',struct('Stack','OFF'),'BkgCM','ON');
    M.Fit;
    E0(9) = M.FitResult.err(2);
   
    %% With Stacked Systematics
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON'),'BkgCM','OFF');
    M.Fit;
    E0Multi(3) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON','FSD','ON'),'BkgCM','OFF');
    M.Fit;
    E0Multi(4) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON','FSD','ON','RF_BF','ON'),'BkgCM','OFF');
    M.Fit;
    E0Multi(5) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON','FSD','ON','RF_BF','ON','RF_EL','ON'),'BkgCM','OFF');
    M.Fit;
    E0Multi(6) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON','FSD','ON','RF_BF','ON','RF_EL','ON','RF_RX','ON'),'BkgCM','OFF');
    M.Fit;
    E0Multi(7) = M.FitResult.err(2);
    
    M.ComputeCM('SysEffect',struct('TCoff_Rad','ON','TCoff_OTHER','ON','TASR','ON','FSD','ON','RF_BF','ON','RF_EL','ON','RF_RX','ON',...
                            'FPDeff','ON'),'BkgCM','OFF');
    M.Fit;
    E0Multi(8) = M.FitResult.err(2);
    
    M.ComputeCM('BkgCM','ON');
    M.Fit;
    E0Multi(9) = M.FitResult.err(2);
    %% for stacked plot
    E0Stacked = zeros(9,1);
    E0Single  = zeros(9,1);
    
    E0Stacked(1) = E0Multi(1);
    E0Single(1)  = E0(1);
    
    for i=2:numel(E0Multi)
        E0Stacked(i) =  E0Multi(i)-E0Multi(i-1);
        E0Single(i)  =   E0(i)-E0(1);
    end
    StackEmpty = zeros(9,1);
    save(savename,'E0','E0Multi','E0Stacked','E0Single','StackEmpty');
end
%% Displayf55 = figure('Name','MultiBarPlot','Renderer','opengl');
f55 = figure('Renderer','opengl');
set(f55, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
b = barh([10;20],[E0Stacked , StackEmpty]','stacked');
b(1).LineStyle = '--';
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';b(5).LineStyle = 'none';b(6).LineStyle = 'none';
b(7).LineStyle = 'none';b(8).LineStyle = 'none';b(9).LineStyle = 'none';
leg_str = {' Statistical';'+ Tritium Acitivity Fluctuation';'+ Final State Distribution';...
    '+ Magnetic Fields';'+ Energy Loss Function'; '+ Column Density + Inel. Cross Section';'+ FPD Efficiency';'+ Background Slope'}; %;'+ Response Function'};
b(1).BarWidth = b(1).BarWidth/2;
%% Single Sys
SingleSysBarWidth = b(1).BarWidth/(2*3);
hold on;
pnone = plot(NaN*[1 1],NaN*[1 1],'Color','w');
%bTC  = barh([7.7;21], [mNuSingle(1),mNuSingle(2); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
%bTC(2).FaceColor =rgb('FireBrick'); bTC(2).LineStyle = 'none'; bTC(2).FaceAlpha = 0.5;
%bTC(1).FaceColor = 'w';%bFSD(1).FaceColor = 'w'; bRF(1).FaceColor = 'w';

bTASR  = barh([7.55;21], [E0Single(1),E0Single(3); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bTASR(2).FaceColor =rgb('FireBrick'); bTASR(2).LineStyle = 'none'; bTASR(2).FaceAlpha = 0.5;
bTASR(1).FaceColor = 'w'; bTASR(1).LineStyle = 'none';

bFSD  = barh([6.62;21], [E0Single(1),E0Single(4); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bFSD(2).FaceColor =rgb('GoldenRod'); bFSD(2).LineStyle = 'none'; bFSD(2).FaceAlpha = 0.5;
bFSD(1).FaceColor = 'w'; bFSD(1).LineStyle = 'none';

bRF_BF  = barh([5.63;21], [E0Single(1),E0Single(5); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bRF_BF(2).FaceColor =rgb('CadetBlue'); bRF_BF(2).LineStyle = 'none'; bRF_BF(2).FaceAlpha = 0.5;
bRF_BF(1).FaceColor = 'w'; bRF_BF(1).LineStyle = 'none';

bRF_EL  = barh([4.57;21], [E0Single(1),E0Single(6); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bRF_EL(2).FaceColor =rgb('PowderBlue'); bRF_EL(2).LineStyle = 'none'; bRF_EL(2).FaceAlpha = 0.9;
bRF_EL(1).FaceColor = 'w'; bRF_EL(1).LineStyle = 'none';

bRF_RX  = barh([3.44;21], [E0Single(1),E0Single(7); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bRF_RX(2).FaceColor =rgb('DarkSlateGray'); bRF_RX(2).LineStyle = 'none'; bRF_RX(2).FaceAlpha = 0.5;
bRF_RX(1).FaceColor = 'w'; bRF_RX(1).LineStyle = 'none';

bFPD  = barh([2.23;21], [E0Single(1),E0Single(8); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bFPD(2).FaceColor =rgb('YellowGreen'); bFPD(2).LineStyle = 'none'; bFPD(2).FaceAlpha = 0.5;
bFPD(1).FaceColor = 'w'; bFPD(1).LineStyle = 'none';

bBkg  = barh([0.935;21], [E0Single(1),E0Single(9); StackEmpty(1:2)'],'stacked','BarWidth',SingleSysBarWidth);
bBkg(2).FaceColor =rgb('Orange'); bBkg(2).LineStyle = 'none'; bBkg(2).FaceAlpha = 0.5;
bBkg(1).FaceColor = 'w'; bBkg(1).LineStyle = 'none';
%% Plot again for nicer statistical border --
b = barh([10;20],[E0Stacked , StackEmpty]','stacked');
b(1).LineStyle = '--';
b(2).LineStyle = 'none';b(3).LineStyle = 'none'; b(4).LineStyle = 'none';b(5).LineStyle = 'none';
b(6).LineStyle = 'none';b(7).LineStyle = 'none';b(8).LineStyle = 'none'; b(9).LineStyle = 'none';
b(1).FaceColor = rgb('White');
b(3).FaceColor = rgb('FireBrick');
b(4).FaceColor = rgb('GoldenRod');
b(5).FaceColor = rgb('CadetBlue');
b(6).FaceColor = rgb('PowderBlue');
b(7).FaceColor = rgb('DarkSlateGray');
b(8).FaceColor = rgb('YellowGreen');
b(9).FaceColor = rgb('Orange');
b(1).BarWidth = b(1).BarWidth/2;

leg_str = {leg_str{:},'Single Contributions','Stat + Tritium Acitivity Fluctuation', 'Stat + Final State Distribution',...
    'Stat + Magnetic Fields','Stat + Energy Loss Function','Stat + Column Density + Inel. Cross Section', 'Stat + FPD Efficiency','Stat + Background Slope'};
leg = legend([b(1),b(3),b(4),b(5),b(6),b(7),b(8),b(9),pnone,bTASR(2),bFSD(2),bRF_BF(2),bRF_EL(2),bRF_RX(2),bFPD(2),bBkg(2)],leg_str{:});
leg.Location = 'north';
leg.NumColumns = 2; legend boxoff
PrettyFigureFormat;
ylim([0 18]);
xlim([0.12,0.26]);
grid on;
set(gca,'FontSize',18);
switch Mode
    case 'Sim'
        xlabel(sprintf('endpoint sensitivity 68.3%% C.L. (eV)'));
    case 'Data'
        xlabel(sprintf('endpoint uncertainty 68.3%% C.L. (eV)'));
end
yticks = [];yticklabels('');

%% Save
plot_dir = strrep(mydir,'results','plots');
if ~exist(plot_dir,'dir')
    system(['mkdir ',plot_dir]);
end
plot_name = [plot_dir,sprintf('FT_E0SensitivityBreakdown_%s.png',Mode)];
print(f55,plot_name,'-dpng','-r450');







