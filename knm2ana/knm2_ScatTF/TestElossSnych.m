RunAnaArg = {'RunNr',RunNr,...
    'fixPar',freePar,...     % free Parameter !!
    'DataType',DataType,...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1,...
    'MosCorrFlag','OFF',...
    'TwinBias_Q',18573.7,...
    'ROIFlag','Default'};

SR = RunAnalysis(RunAnaArg{:});

qu = 18575;%SR.ModelObj.qU(1);
E  = 18575;
Theta = linspace(0,0.88,100);

% new: effective model
CommonArg = {'qU',qu,'E',E,'Theta_rad',Theta,'Bs_T',SR.ModelObj.WGTS_B_T,'Bt_T',3.6};
ElossPara = SR.ModelObj.GetElossSync(CommonArg{:},'Mode','Combi');
ElossPara_s = SR.ModelObj.GetElossSync(CommonArg{:},'Mode','Source'); % source only
ElossPara_t = SR.ModelObj.GetElossSync(CommonArg{:},'Mode','Transp'); % transport only

% OLD: Samak Simulation Providing Synchrotron Loss in WGTS
PinchAngle =   [0.017453    0.034907     0.05236    0.069813    0.087266     0.10472     0.12217     0.13963     0.15708     0.17453     0.19199     0.20944     0.22689     0.24435      0.2618     0.27925     0.29671     0.31416     0.33161     0.34907     0.36652     0.38397     0.40143     0.41888     0.43633     0.45379     0.47124     0.48869     0.50615      0.5236     0.54105     0.55851     0.57596     0.59341     0.61087     0.62832     0.64577     0.66323     0.68068     0.69813     0.71558     0.73304     0.75049     0.76794      0.7854     0.80285      0.8203     0.83776     0.85521     0.87266];
SeLoss     =   [8.4048e-06   3.363e-05  7.5709e-05   0.0001347  0.00021067  0.00030374  0.00041402  0.00054166  0.00068685  0.00084979   0.0010307   0.0012299   0.0014476   0.0016842   0.0019401   0.0022157   0.0025114   0.0028278   0.0031655   0.0035251   0.0039073   0.0043129   0.0047428    0.005198   0.0056796   0.0061889   0.0067271   0.0072959   0.0078971   0.0085326   0.0092046   0.0099157    0.010669    0.011468    0.012316    0.013218    0.014179    0.015206    0.016306    0.017488    0.018763    0.020145     0.02165    0.023299    0.025119    0.027147    0.029428    0.032031    0.035047    0.038618];
DEsyn   = @(angleRad) interp1(PinchAngle,SeLoss,angleRad,'spline','extrap'); % data go till 50 deg, but computation needs 51 deg or so --> spline extrapolation
ElossOld = DEsyn(Theta);

%% plot comparison
figure('Units','normalized','Position',[0.1,0.1,0.4,0.4]);
p1 =plot(Theta,ElossOld,'-','LineWidth',2,'Color',rgb('FireBrick'));
hold on;
% pData = errorbar(PinchAngle,SeLoss,zeros(numel(SeLoss),1),'o',...
%     'LineWidth',1,'Color',rgb('FireBrick'),...
%     'MarkerFaceColor',rgb('FireBrick'),...
%     'CapSize',0);
p2 = plot(Theta,ElossPara,'-','LineWidth',2,'Color',rgb('DodgerBlue'));
p2_s = plot(Theta,ElossPara_s,':','LineWidth',2,'Color',rgb('LimeGreen'));
p2_t = plot(Theta,ElossPara_t,'-.','LineWidth',2,'Color',rgb('PowderBlue'));
PrettyFigureFormat;
hold off
xlabel(sprintf('\\theta (rad)'));
ylabel('Synchrotron energy loss (eV)');
leg = legend([p1,p2,p2_t,p2_s],'Interpolation (old)','Parametrization (new)',...
    'Parametrization transport (new)','Parametrization source (new)',...
    'EdgeColor',rgb('Silver'),'Location','northwest');
xlim([0 max(Theta)])
