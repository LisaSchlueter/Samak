% KNM2 Background Fit
% Extract the Slope
% Uniform

DiagnosticFlag = 'OFF';

% Read KNM2 data
range         = 5;
RunAnaArg = {...
    'RunList','KNM2_Prompt',...        % all KNM2 golden runs
    'DataType','Real',...
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'ROIFlag','Default',...
    };%
NonPoissonScaleFactor = 1.112;
A                     = MultiRunAnalysis(RunAnaArg{:});
A.exclDataStart       = A.GetexclDataStart(range);

%% Fit Slope
nbin  = numel(A.RunData.qU(A.exclDataStart:end));
hv    = A.RunData.qU(A.exclDataStart:end);
c     = A.RunData.TBDIS(A.exclDataStart:end)./A.RunData.qUfrac(A.exclDataStart:end)/A.RunData.TimeSec;
ce    = NonPoissonScaleFactor*A.RunData.TBDISE(A.exclDataStart:end)./A.RunData.qUfrac(A.exclDataStart:end)/A.RunData.TimeSec;
% column1: hv in kV
% column2: rate in mcps
% column3: rate error in mcps 
Data = [hv/1e3 c*1e3 ce*1e3];
%Data = Data([1 7],:);
% Linear Model
parnames = {'baseline','slope'}; 
ParIni   = [mean(Data(:,2)) 0]; 

% Case 1) Fit Overall Range
% Fit
Args = {ParIni, Data, '-c', 'min; minos'};
[ par, err, chi2min, errmat ] = fminuit('Chi2Gauss',Args{:});
fprintf('=========================== Fit results ========================\n');
fprintf('  baseline = %.3f ± %.3f\n',par(1),err(1));
fprintf('  slope    = %.3f ± %.3f\n',par(2),err(2));
fprintf('  chi2min  = %.2f for ndof = %d\n',chi2min,nbin-numel(ParIni));
fprintf('================================================================\n');

% Figure 1
close all
figure(1)
set(gcf, 'Position',  [100, 100, 1000, 500])
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),'ks',...
       'MarkerSize',8,'MarkerFaceColor',.8*[1 1 1],'LineWidth',2);
hold on;
hfit    = plot(Data(:,1),Model(par,Data(:,1)),'LineWidth',3);
%hfiteu =  plot(Data(:,1),Model(par-[0 err(2)],Data(:,1)),'LineWidth',1,'LineStyle','--');
%hfited =  plot(Data(:,1),Model(par+[0 err(2)],Data(:,1)),'LineWidth',1,'LineStyle','--');
hold off;
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('rate (mcps)','FontSize',14) 
PrettyFigureFormat;
legend([hdata hfit],...
    sprintf('KNM2 Data: %.0f counts in %.0f sec',sum(A.RunData.TBDIS(A.exclDataStart:end)),sum(A.RunData.qUfrac(A.exclDataStart:end).*A.RunData.TimeSec)),...
    sprintf('Fit \n baseline=%0.3f \\pm %0.3f mcps \n slope=%0.3f \\pm %0.3f mcps/keV \n relative slope = %0.2f \\pm %0.2f %%/keV',...
    par(1),err(1),par(2),err(2),par(2)/par(1)*100,err(2)/par(1)*100),'Location','NorthEast')
grid on
publish_figurePDF(1,'BKG_KNM2.pdf');

switch DiagnosticFlag
    case 'OFF'
        return;
end

%% Fit Diagnostics
% Put the variables into a table, naming them appropriately
tbl = table(Data(:,1)-18.5737,Data(:,2),Data(:,3),'VariableNames',{'qUKEV','RateMCPS','RateErrorMCPS'});
% Specify and carry out the fit
mdl = fitlm(tbl,'RateMCPS ~ 1 + qUKEV','Weights',1./((Data(:,3)).^2),'RobustOpts','off','Exclude',[]);

% Leverage
figure(300)
plotDiagnostics(mdl,'leverage','Color','b','MarkerSize',10,'Marker','s','LineWidth',2);
legend('show') % Show the legend
PrettyFigureFormat

% Cook Distance
figure(301)
plotDiagnostics(mdl,'cookd','Color','b','MarkerSize',10,'Marker','s','LineWidth',2);
legend('show') % Show the legend
PrettyFigureFormat

% Dfbetas
figure(302)
plotDiagnostics(mdl,'Dfbetas','Color','b','MarkerSize',10,'Marker','s','LineWidth',2);
legend('show') % Show the legend
PrettyFigureFormat

%% Change of Slope
figure(555)
set(gcf, 'Position',  [100, 100, 1000, 700])
subplot(2,1,1)
hdata = errorbar(Data(:,1)-18.5737,Data(:,2),Data(:,3),'ks',...
    'MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('Rate (mcps)','FontSize',14) 
grid on
PrettyFigureFormat
subplot(2,1,2)
h=plot(Data(:,1)-18.5737,mdl.Diagnostics.Dfbetas(:,2) .* mdl.Coefficients.SE(2) ./ mdl.Coefficients.Estimate(2)*100,'ks',...
    'MarkerSize',8,'MarkerFaceColor',rgb('IndianRed'),'LineWidth',2);
xlabel('Retarding Potential (keV)','FontSize',14) 
ylabel('\Delta Slope (%)','FontSize',14) 
grid on
legend([h],sprintf('all data slope \n %0.3f \\pm %0.3f mcps/keV',par(2),err(2)),'Location','SouthEast')
PrettyFigureFormat

