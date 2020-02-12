function [Q_i_Fit, chi2min, dof, E0min, E090low, E090up] = E0Scan_SensitivityStudy_FTKATRIN(varargin)
addpath(genpath('../../../Samak2.0'));
% ------------------------------------------------------------------------------------------------
% Goal: 
% Obtain sensitivity on effective endpoint 
% Method: Scan
% Simulate KATRIN 1 time (asimov): First Tritium (June 2018) settings
% Options: SysEffect, energy range, plot,...
% Fit with statistic only or with a covariance matrix (E0,N,B free) and a fixed E0 (scan over E0)
% Plot chi^2 as a function of E0
% Compute sensitivity on E0 (chi2(min)+2.7)
% ------------------------------------------------------------------------------------------------
% Lisa Schlüter, September 18 (TUM/MPP)
%---------------------------------------parser----------------------------------------------------------%
p = inputParser;
p.addParameter('chi2','chi2CM',@(x)ismember(x,{'chi2Stat','chi2CM'}));
p.addParameter('SysEffect','all',@(x)ischar(x)); %1 sys effect e.g. FSD
p.addParameter('exclDataStart',7,@(x)isfloat(x));
p.addParameter('Q_i_Fit',-0.25:0.005:0.25,@(x)all(isfloat(x)));
p.addParameter('plotFit','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('RunList','StackCD100all',@(x)ischar(x)); %

p.parse(varargin{:});

chi2          = p.Results.chi2;
SysEffect     = p.Results.SysEffect;
exclDataStart = p.Results.exclDataStart;
Q_i_Fit       = p.Results.Q_i_Fit;
plotFit       = p.Results.plotFit;
RunList       = p.Results.RunList;
%---------------------------------------parser end----------------------------------------------------------%
if strcmp(chi2,'chi2Stat') %for labeling 
    SysEffect = '';       
end
nQi = numel(Q_i_Fit);

%Init Model
MRA = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat',...%chi2 set later
                        'fixPar','1 2 5 6','exclDataStart',exclDataStart);
belowE0 = MRA.ModelObj.Q_i-MRA.ModelObj.qU(exclDataStart)-1; % for labeling
MRA.InitModelObj_Norm_BKG;              % Set Background to level of real data
MRA.RunData.TBDIS = MRA.ModelObj.TBDIS; % Give Simulated Data to MRA


%Covariance Matrix
MRA.chi2 = chi2; 
if ~strcmp(chi2,'chi2Stat')
    if strcmp(SysEffect,'RF')
        MRA.ComputeCM('SysEffects',struct('RF_EL','ON','RF_BF','ON','RF_RX','ON'),'StackCM','OFF',...
            'DataDriven','OFF','InitNormFit','ON');
    elseif strcmp(SysEffect,'TC')
        MRA.ComputeCM('SysEffects',struct('TCoff_OTHER','ON','TCoff_RAD','ON'),'StackCM','OFF',...
            'DataDriven','OFF','InitNormFit','ON');
    elseif strcmp(SysEffect,'all')
        MRA.ComputeCM('StackCM','ON','DataDriven','OFF','InitNormFit','ON');
    else
        MRA.ComputeCM('SysEffects',struct(SysEffect,'ON'),'StackCM','OFF',...
            'DataDriven','OFF','InitNormFit','ON');
    end
else
    MRA.InitModelObj_Norm_BKG;          % Init to Simulation
end

FitResult = cell(nQi,1);
chi2min   = zeros(nQi,1);

%Start Scan
progressbar(sprintf('Endpoint Scan %s- FT Sensitivity Study',SysEffect));
for i=1:nQi
    progressbar(i/nQi);
   % MRA.ModelObj.Q_i = 18573.7+Q_i_Fit(i);  
    MRA.Fit;
    FitResult{i} = MRA.FitResult;
    chi2min(i) = MRA.FitResult.chi2min;
end

%Compute Sensitivity (90% C.L.)
E090low =interp1(chi2min(1:ceil(nQi/2)),Q_i_Fit(1:ceil(nQi/2)),min(chi2min)+2.7,'spline');
E090up = interp1(chi2min(ceil(nQi/2):end),Q_i_Fit(ceil(nQi/2):end),min(chi2min)+2.7,'spline');
E0min = interp1(chi2min,Q_i_Fit,min(chi2min),'spline');

% Check if E0(90% .C.L.) is inside scanning range. Redo scan if not;
if E090low<=min(Q_i_Fit) || E090up>=max(Q_i_Fit)
    Q_i_Lim =round(max([E090low E090up]),1);
    Q_i_Fit= -Q_i_Lim:0.01:Q_i_Lim;
    nQi = numel(Q_i_Fit);
    if nQi>=200
        fprintf(2,'Redo Scan with different Q_i_Fit. 90%%C.L. not inside scan range \n');
        return;
    end
    progressbar('REDO Endpoint Scan - FT Sensitivity Study');
    for i=1:nQi
        progressbar(i/nQi);
        MRA.ModelObj.Q_i = 18573.7+Q_i_Fit(i);
        MRA.Fit;
        FitResult{i} = MRA.FitResult;
        chi2min(i) = MRA.FitResult.chi2min;
    end
    E090low =interp1(chi2min(1:ceil(nQi/2)),Q_i_Fit(1:ceil(nQi/2)),min(chi2min)+2.7,'spline');
    E090up = interp1(chi2min(ceil(nQi/2):end),Q_i_Fit(ceil(nQi/2):end),min(chi2min)+2.7,'spline');
    E0min = interp1(chi2min,Q_i_Fit,min(chi2min),'spline');
end
dof = FitResult{1}.dof;
%% plot
if strcmp(plotFit,'ON')
f3 = figure(3);
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
x = [E090low, E090up];
y = [min(chi2min),min(chi2min)+2.7];
pScan = plot(Q_i_Fit,chi2min,'LineWidth',3,'Color',rgb('CadetBlue'));
hold on;
pLine = plot(x,(min(chi2min)+2.7).*ones(numel(x),1),'--','LineWidth',3,'Color',rgb('IndianRed'));
plot(E090low.*ones(numel(y),1),y,'--','LineWidth',3,'Color',rgb('IndianRed'));
plot(E090up.*ones(numel(y),1),y,'--','LineWidth',3,'Color',rgb('IndianRed'));
xlabel('\Delta E0_{eff} (eV)')
ylabel(['\chi^2 (',num2str(dof),'dof)']);
xlim([min(Q_i_Fit) max(Q_i_Fit)]);
ylim([0 max(chi2min)]);
scan_leg = ['\Delta E0_{eff} = ',sprintf('%.3f - %.3f + %.3f eV (90%% C.L.)\n',E0min,abs(E090low), E090up)];
line_leg = sprintf('\\chi2_{min}+2.7');
legend([pScan, pLine],scan_leg,line_leg,'Location','northwest'); legend boxoff;
PrettyFigureFormat;
set(gca,'FontSize',18);
title(sprintf('Samak %s %s Fit to Simulation (asimov) \n - KATRIN First Tritium %.1f days - MTD: %s - BKG %.0fmcps - %.0feV range',...
    chi2,SysEffect,MRA.ModelObj.TimeSec/(60*60*24),MRA.StackFileName,MRA.ModelObj.BKG_RateSec*1e3,belowE0));
%% save
if ~exist('./plots/png','dir')
    mkdir ./plots/png/
    mkdir ./plots/pdf/
    mkdir ./plots/fig/
end
save_name = sprintf('SensitivityFT_E0Scan_%s%s_%s_%.0feV',chi2,SysEffect,MRA.StackFileName,belowE0);
print(f3,['./plots/png/',save_name,'.png'],'-dpng');
savefig(f3,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(f3,['./plots/pdf/',save_name,'.pdf']);
end
end