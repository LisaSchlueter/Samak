% Compute neutrino mass sensitivity before data taking with fake MC run

% Simulation Parameters
Flag  = 'KNM1';
NuMassSquare = 1;    % eV2
DeltaE0      = 0.08; % eV
range        = [-50 -1];

% init files
Debug = 'OFF';
PlotStyle = { 'o',... 
    'MarkerSize',16,...
    'MarkerFaceColor',...
    rgb('SkyBlue'),... 
    'LineWidth',2,...
    'Color',rgb('IndianRed')};
LocalFontSize = 24;
                    
% Common Arguments
CommonArg = {'RunNr',1,...% has no meaning
    'DataType','Fake',...
    'FSDFlag','BlindingKNM2',...
    'fixPar','Norm E0 BBkg',...
    'ELossFlag','KatrinT2',...
    'AnaFlag','StackPixel',...
    'chi2','chi2Stat',...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    };


% Start
switch Flag
    case 'KNM1'
        InitFile = @ref_FakeRun_KNM1_CD22_23days;
        SpecArg = { 'exclDataStart',14,... % 40eV range
            };
    case 'KNM2'
        InitFile = @ref_FakeRun_KNM2_CD84_30days;
        SpecArg = { 'exclDataStart',12,... % 40eV range
            };
end


%% Zero Neutrino Mass
R0 = RunAnalysis(CommonArg{:},SpecArg{:},'FakeInitFile',InitFile);
R0.ModelObj.mnuSq_i=0;

%% 1 eV2 Neutrino Mass Squared
R1 = RunAnalysis(CommonArg{:},SpecArg{:},'FakeInitFile',InitFile);
R1.ModelObj.mnuSq_i=NuMassSquare;

%% 0 eV2 Neutrino for convolution
R2 = RunAnalysis(CommonArg{:},SpecArg{:},'FakeInitFile',InitFile);
R2.ModelObj.mnuSq_i=0;

%% 0 eV2 Neutrino with shifted endpoint
R3 = RunAnalysis(CommonArg{:},SpecArg{:},'FakeInitFile',InitFile);
R3.ModelObj.mnuSq_i=0;

%% Background
switch Flag
    case 'KNM1'
        R3.ModelObj.BKG_RateSec_i = 0.29; % cps
        R2.ModelObj.BKG_RateSec_i = 0.29; % cps
        R1.ModelObj.BKG_RateSec_i = 0.29; % cps
        R0.ModelObj.BKG_RateSec_i = 0.29; % cps
        m2OnesigmaSensitivity     = 1;  % eV2
        mROI                      = 12;   % eV
    case 'KNM2'
        R3.ModelObj.BKG_RateSec_i = 0.21; % cps
        R2.ModelObj.BKG_RateSec_i = 0.21; % cps
        R1.ModelObj.BKG_RateSec_i = 0.21; % cps
        R0.ModelObj.BKG_RateSec_i = 0.21; % cps
        m2OnesigmaSensitivity     = 0.5;  % eV2
        mROI                      = 9; % eV
end

%% Compute Spectra
R0.ModelObj.ComputeTBDDS; R0.ModelObj.ComputeTBDIS;
R1.ModelObj.ComputeTBDDS; R1.ModelObj.ComputeTBDIS;
R2.ModelObj.ComputeTBDDS; R2.ModelObj.ComputeTBDIS;

%% R2: Spectral Convolution
sigmaBroad=sqrt(abs(R1.ModelObj.mnuSq_i))/sqrt(2);
cd=gaussfilt(R0.ModelObj.Te,R0.ModelObj.TBDDS,sigmaBroad);
R2.ModelObj.TBDDS=cd;R2.ModelObj.ComputeTBDIS;

%% R3: Shifted Endpoint
%DeltaE0 = m2OnesigmaSensitivity ./ 2 ./ mROI; % eV
R3.ModelObj.Q_i = R3.ModelObj.Q_i  + DeltaE0; %eV
R3.ModelObj.ComputeTBDDS; R3.ModelObj.ComputeTBDIS;

switch Debug
    case 'ON'
%% Plot Differential Spectrum
figure0 = figure();
eme0 = R0.ModelObj.Te-R0.ModelObj.Q_i+1.72;
set(figure0, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
s1=subplot(2,1,1)
h0=plot(eme0,R0.ModelObj.TBDDS ,'Color','Black','LineWidth',2);
hold on
h1=plot(eme0,R1.ModelObj.TBDDS,'Color','Red','LineWidth',2);
hold off
legend([h0, h1],...
    sprintf('%s : m=%.1f eV^2',Flag,R0.ModelObj.mnuSq),...
    sprintf('%s : m=%.1f eV^2',Flag,R1.ModelObj.mnuSq));
  legend boxoff;
  set(gca, 'YScale', 'lin');
xlabel(sprintf('E  -  %.1f  (eV)',R0.ModelObj.Q_i),'FontSize',20);
ylabel('Differential Spectrum (cps)','FontSize',20);
PrettyFigureFormat('FontSize',LocalFontSize);
s2=subplot(2,1,2)
h0=plot(eme0,R1.ModelObj.TBDDS./R0.ModelObj.TBDDS ,'Color',rgb('IndianRed'),'LineWidth',2);
xlabel(sprintf('E  -  %.1f  (eV)',R0.ModelObj.Q_i),'FontSize',20);
ylabel('Ratio','FontSize',20);
            PrettyFigureFormat('FontSize',LocalFontSize); 
linkaxes([s1,s2],'x');
xlim(range);
end

switch Debug
    case 'ON'
%% Plot Integral Spectrum
% ./(-R0.ModelObj.qU+R0.ModelObj.Q_i).^3
figure1 = figure();
set(figure1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);
h0=plot(R0.ModelObj.qU-R0.ModelObj.Q_i,...
    R0.ModelObj.TBDIS./R0.ModelObj.TimeSec./R0.ModelObj.qUfrac,...
    'Color','Black','LineWidth',2);
hold on
b0=line([-100 100],[R0.ModelObj.BKG_RateSec_i R0.ModelObj.BKG_RateSec_i],'Color',rgb('IndianRed'),'LineStyle','--');
s0=plot(R0.ModelObj.qU-R0.ModelObj.Q_i,...
    (R0.ModelObj.TBDIS-R0.ModelObj.BKG_RateSec_i.*R0.ModelObj.TimeSec.*R0.ModelObj.qUfrac)./R0.ModelObj.TimeSec./R0.ModelObj.qUfrac,...
    'LineWidth',2,'Color',rgb('DarkGreen'),'LineStyle','--');
h1=plot(R0.ModelObj.qU-R0.ModelObj.Q_i,...
    R1.ModelObj.TBDIS./R0.ModelObj.TimeSec./R0.ModelObj.qUfrac,...
    'Color','Red','LineWidth',2);
hold off
legend([h0, h1],...
    sprintf('%s : m=%.1f eV^2',Flag,R0.ModelObj.mnuSq),...
    sprintf('%s : m=%.1f eV^2',Flag,R1.ModelObj.mnuSq));
  legend boxoff;
  set(gca, 'YScale', 'log');
xlim([-40 +5]);ylim([0.1 1000]);
xlabel(sprintf('mean retarding potential < qU >  -  %.1f  (eV)',R0.ModelObj.Q_i),'FontSize',20);
ylabel('Integral Spectrum (cps)','FontSize',20);
            PrettyFigureFormat('FontSize',LocalFontSize); 
end

%% Plot Integral Spectrum Ratio to Zero Neutrino Mass
figure2 = figure();
set(figure2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);
ratio0 = line([-100 100],[1,1],'Color','Black','LineStyle','--');
hold on
ratio1=errorbar(R0.ModelObj.qU-R0.ModelObj.Q_i,...
    R1.ModelObj.TBDIS./R0.ModelObj.TBDIS,...
    sqrt(R1.ModelObj.TBDIS)./R0.ModelObj.TBDIS,...
    PlotStyle{:});
hold off
xlabel(sprintf('mean retarding potential < qU >  -  %.1f  (eV)',R0.ModelObj.Q_i),'FontSize',20);
ylabel('Ratio','FontSize',20);
legend([ratio1],...
    sprintf('%s : m=%.1f eV^2',Flag,R1.ModelObj.mnuSq),'Location','NorthWest');
  legend boxoff;
  xlim(range);
            PrettyFigureFormat('FontSize',LocalFontSize); 

%% Convolution Differential Spectrum
switch Debug
    case 'ON'
        figure(999)
d  = R0.ModelObj.TBDDS;
s1=subplot(2,1,1)
plot(R0.ModelObj.Te,d,'Color',rgb('IndianRed'),'LineStyle','--');
hold on
plot(R0.ModelObj.Te,cd,'Color',rgb('DarkGreen'),'LineStyle','-');
hold off
set(gca, 'YScale', 'log');
            PrettyFigureFormat('FontSize',LocalFontSize); 
s2=subplot(2,1,2)
plot(R0.ModelObj.Te,cd./d,'Color',rgb('DarkGreen'),'LineStyle','-');
set(gca, 'YScale', 'log');
            PrettyFigureFormat('FontSize',LocalFontSize); 
linkaxes([s1,s2],'x');
end

%% Convolution Integram Spectrum
switch Debug
    case 'ON'
        figure(9999)
s1=subplot(2,1,1)
plot(R0.ModelObj.qU-R0.ModelObj.Q_i,R0.ModelObj.TBDIS./R0.ModelObj.TimeSec./R0.ModelObj.qUfrac,'Color',rgb('IndianRed'),'LineStyle','--');
hold on
plot(R0.ModelObj.qU-R0.ModelObj.Q_i,R2.ModelObj.TBDIS./R0.ModelObj.TimeSec./R0.ModelObj.qUfrac,'Color',rgb('DarkGreen'),'LineStyle','-');
hold off
set(gca, 'YScale', 'log');
            PrettyFigureFormat('FontSize',LocalFontSize); 
s2=subplot(2,1,2)
plot(R0.ModelObj.qU-R0.ModelObj.Q_i,R2.ModelObj.TBDIS./R0.ModelObj.TBDIS,'Color',rgb('DarkGreen'),'LineStyle','-');
set(gca, 'YScale', 'log');
            PrettyFigureFormat('FontSize',LocalFontSize); 
linkaxes([s1,s2],'x');
end

%% Plot Integral Spectrum Ratio to Zero Neutrino Mass
figure3 = figure();
set(figure3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);
ratio0 = line([-100 100],[1,1],'Color','Black','LineStyle','--');
hold on
ratio1=errorbar(R0.ModelObj.qU-R0.ModelObj.Q_i,...
    R2.ModelObj.TBDIS./R0.ModelObj.TBDIS,...
    sqrt(R2.ModelObj.TBDIS)./R0.ModelObj.TBDIS,...
    PlotStyle{:});
hold off
xlabel(sprintf('mean retarding potential < qU >  -  %.1f  (eV)',R0.ModelObj.Q_i),'FontSize',20);
ylabel('Ratio','FontSize',20);
legend([ratio1],...
    sprintf('%s : Unknown Gaussian Broadening \\sigma = %.1f eV',Flag,sigmaBroad),'Location','NorthWest');
  legend boxoff;
  xlim(range);
            PrettyFigureFormat('FontSize',LocalFontSize); 

%% Plot Integral Spectrum Ratio to Zero Neutrino Mass
figure4 = figure();
set(figure4, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);
ratio0 = line([-100 100],[1,1],'Color','Black','LineStyle','--');
hold on
ratio1=errorbar(R0.ModelObj.qU-R0.ModelObj.Q_i,...
    R3.ModelObj.TBDIS./R0.ModelObj.TBDIS,...
    sqrt(R3.ModelObj.TBDIS)./R0.ModelObj.TBDIS,...
    PlotStyle{:});
hold off
xlabel(sprintf('mean retarding potential < qU >  -  %.1f  (eV)',R0.ModelObj.Q_i),'FontSize',20);
ylabel('Ratio','FontSize',20);
legend([ratio1],...
    sprintf('%s : Endpoint Shift w.r.t true value \\delta E_0 = %.2f eV',Flag,DeltaE0),'Location','NorthWest');
  legend boxoff;
  xlim(range);
            PrettyFigureFormat('FontSize',LocalFontSize); 
