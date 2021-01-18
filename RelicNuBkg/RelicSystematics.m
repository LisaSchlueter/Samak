startup;
Recompute = 'OFF';
B=RelicNuDebug('Params','KNM1');
matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Misc/')];
%B.SystBreakdown('TwinBias_mnuSq',0,'Plot','OFF','Recompute',Recompute);
B.SystBreakdown('TwinBias_mnuSq',1,'Plot','OFF','Recompute',Recompute);
B.SystBreakdown('TwinBias_mnuSq',5,'Plot','OFF','Recompute',Recompute);
load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq0_range40.mat']);
Y0=Y;
load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq1_range40.mat']);
Y1=Y;
load([matFilePath,'SensitivityBreakdown_KNM1_mnuSq5_range40.mat']);
Y5=Y;

bar(X,[Y0;Y1;Y5]);
ylabel('\eta');
legend('m_{\nu}^{2}=0 eV^{2}','m_{\nu}^{2}=1 eV^{2}','m_{\nu}^{2}=5 eV^{2}');
legend boxoff;
PrettyFigureFormat;