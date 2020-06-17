%% Sterile theoretical spectra

E0 = 18570;          % Endpoint in eV

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
        
times = (R.ModelObj.qUfrac*R.ModelObj.TimeSec);

R.ModelObj.ComputeTBDDS();
YD=R.ModelObj.TBDDS;
R.ModelObj.ComputeTBDIS();

YI = R.ModelObj.TBDIS;
YI = YI./times;

qU=R.ModelObj.qU;   % Energy axis
qU=qU-E0;

YI=YI(qU>-90);
qUc=qU(qU>-90);
%% Sterile
%sterile_mass = 5;     % Sterile neutrino mass in eV
mixing_angle = 0.1;    % sin(th4)

for sterile_mass=[0.5 2 5 15 30]
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
%   Rs.ModelObj.mnuSq_i  = 0.9;

    Rs.ModelObj.ComputeTBDDS();
    YDs = Rs.ModelObj.TBDDS;
    Rs.ModelObj.ComputeTBDIS();

    YIs = Rs.ModelObj.TBDIS;
    err = sqrt(YIs)./times;
    YIs = YIs./times;
    
    YIs=YIs(qU>-90);
    err=err(qU>-90);
    
    RSP=YIs./YI;
    errorbar(qUc,RSP,err./YI)
    hold on;
end

%% ====== Plots ======

%% TBDS
%plot(qUc,YI)
% hold on;
% plot(qUc,YIs)

%% Ratio
% RSP=YIs./YI;        % Ratio
% plot(qUc,RSP)

%% Plot parameters
xlabel('qU - E_0 (eV)');
ylabel('Ratio'); %'Rate (cps)'
legend({'0.5 eV','2 eV','5 eV','15 eV','30 eV'},'Location','southwest') %{'5 eV','10 eV','20 eV','30 eV'}{'Without sterile','With sterile'}
PrettyFigureFormat;
set(gca, 'YScale', 'log');