%% Sterile theoretical spectra

E0 = 18570;                                      % Endpoint in eV
sterile_mass = 20;                               % Sterile neutrino mass in eV
mixing_angle = 0.1;                              % sin2(2th4)
%mixing_angle = (1-sqrt(1-mixing_angle))/2;      % sin2(th4)

%% No sterile

R = MultiRunAnalysis('RunList','KNM1',...
            'chi2','chi2Stat',...
            'DataType','Twin',...
            'fixPar','',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag','Sibille0p5eV',...
            'ELossFlag','KatrinT2',...
            'SysBudget',22);

% Global variables
times = (R.ModelObj.qUfrac*R.ModelObj.TimeSec);

qU = R.ModelObj.qU;     % Energy axis
qU = qU-E0;

% Spectrum
R.ModelObj.ComputeTBDDS();
YD=R.ModelObj.TBDDS;
R.ModelObj.ComputeTBDIS();

YI = R.ModelObj.TBDIS;
YI = YI./times;

%% Sterile theoretical

Rs = MultiRunAnalysis('RunList','KNM1',...
            'chi2','chi2Stat',...
            'DataType','Twin',...
            'fixPar','',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag','Sibille0p5eV',...
            'ELossFlag','KatrinT2',...
            'SysBudget',22);
        
Rs.ModelObj.mnu4Sq_i = sterile_mass^2;
Rs.ModelObj.sin2T4_i = mixing_angle;

% Spectrum
Rs.ModelObj.ComputeTBDDS();
YDs=Rs.ModelObj.TBDDS;
Rs.ModelObj.ComputeTBDIS();

IS = Rs.ModelObj.TBDIS;     % Counts
YIs = IS./times;

%% Sterile "data"

YIsd = IS;

% Error
err  = sqrt(YIsd);

% Fluctuations (data sim)
YIsd = YIsd + err.*randn(length(YIsd),1);

YIsd = YIsd./times;

% Error bar
err  = err./times;
err  = err./YI;

%% Constraining to -90eV
YIsd=YIsd(qU>-90);
YIs=YIs(qU>-90);
YI=YI(qU>-90);
err=err(qU>-90);
qUc=qU(qU>-90);

%% Ratio
figure;
RSP  = YIs./YI;
RSPd = YIsd./YI;
plot(qUc,RSP)
hold on;
errorbar(qUc,RSPd,err,'o')

%% Plot parameters
xlabel('qU - E_0 (eV)');
ylabel('Ratio'); %'Rate (cps)'
legend({'Theoretical','Data'},'Location','southwest') %{'5 eV','10 eV','20 eV','30 eV'}{'Without sterile','With sterile'}
PrettyFigureFormat;
set(gca, 'YScale', 'log');