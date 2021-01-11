startup;
Q = RelicNuDebug('Params','KNM1');
eta_Chi2 = zeros(1,21);
for i=1:21
    mNu = (i-1)/10;
    Q.Chi2Twin('Recompute','ON','Plot','OFF','Syst','ON','fitPar','mNu E0 Norm Bkg','DeltaChi2',1,'TwinBias_mnuSq',mNu,'NetaBins',2,'etarange',11,'etafactor',2);
    eta_Chi2(i) = Q.etaSensitivity;
end
[mNu,eta] = Q.EtaFit('Recompute','ON');
save('./Results_mNuFree.mat','eta','eta_Chi2');
C=plot(mNu,eta_Chi2,'LineWidth',2);
hold on;
F=plot(mNu,eta,'LineWidth',2);
xlabel('m_{\nu}^{2}');
ylabel('\eta');
legend('\chi^{2} Scan','\eta Fit');
legend boxoff;
PrettyFigureFormat;