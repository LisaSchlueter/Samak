% B = RelicNuDebug('Params','KNM1');
% B.Chi2Scan_2D('Recompute','ON','Netabins',20,'Nmnubins',20,'etalow',5e10,'etahigh',3.2e11,'mnulow',0.1,'mnuhigh',3,'TwinBias_mnuSq',1);
% B.SystBias;

% A=ref_RelicNuBkg_KNM1('TTFSD','OFF','HTFSD','OFF','DTFSD','OFF','mnuSq_i',1);
% B=ref_RelicNuBkg_KNM1;%('mnuSq_i',1);
% a=semilogy(A.Te-A.Q,A.TBDDS,'LineWidth',2);
% hold on;
% b=semilogy(A.Te-A.Q,A.TBDDS_R,'LineWidth',2);
% c=semilogy(B.Te-B.Q,B.TBDDS,'LineWidth',2);
% d=semilogy(B.Te-B.Q,B.TBDDS_R,'LineWidth',2);
% ylabel('Rate per Energy Bin');
% xlabel('Energy - E_{0} (eV)');
% legend([a b c d],'\beta decay, FSD off','C\nuB, FSD off','\beta decay, FSD on','C\nuB, FSD on');
% legend boxoff;
% PrettyFigureFormat;

% Q   = 18575;
% mnu = 1;
% me  = 511;
% Te  = linspace(Q-10,Q+5,1000);
% B   = ref_RelicNuBkg_KNM1;
% TBDDS1 = (((Q-Te)>=0).*sqrt((me+Te).^2+me^2).*(Te+me).*(Q-Te).*(((Q-Te).^2-mnu^2) >= 0).*((Q-Te).^2-mnu^2).^0.5)+1e9.*exp(-(Te-Q-mnu).^2/(2*0.094^2));
%TBDDS2 = 0.0108.*(((Q-0.3-Te)>=0).*sqrt((me+Te).^2+me^2).*(Te+me).*(Q-0.3-Te).*(((Q-0.3-Te).^2-mnu^2) >= 0).*((Q-0.3-Te).^2-mnu^2).^0.5)+exp(-(Te-Q-0.3-2*mnu)/(2*0.94));
%TBDDS3 = 0.0353.*(((Q-0.8-Te)>=0).*sqrt((me+Te).^2+me^2).*(Te+me).*(Q-0.8-Te).*(((Q-0.8-Te).^2-mnu^2) >= 0).*((Q-0.8-Te).^2-mnu^2).^0.5)+exp(-(Te-Q-0.8-2*mnu)/(2*0.94));
% a = semilogy(Te-Q,0.9526.*TBDDS1,'LineWidth',2);
% hold on;
% b = plot(Te-Q-0.3,0.0108.*TBDDS1,'LineWidth',2);
% c = plot(Te-Q-0.8,0.0353.*TBDDS1,'LineWidth',2);
% xlim = [-4,1.5];
% ylim = [1e6,5.3e9];
%TTexP = B.TTexP;
%TTexE = B.TTexE;
%DTexP = B.DTexP;
%DTexE = B.DTexE;
%HTexP = B.HTexP;
%HTexE = B.HTexE;
%TBDDS2 = zeros(1,numel(TBDDS1));
%for i=1:numel(TTexP)
%    TBDDS2 = TBDDS2 + (TTexP(i).*(((Q-TTexE(i)-Te)>=0).*sqrt((me+Te).^2+me^2).*(Te+me).*(Q-TTexE(i)-Te).*(((Q-TTexE(i)-Te).^2-mnu^2) >= 0).*((Q-TTexE(i)-Te).^2-mnu^2).^0.5));
%end
%b = semilogy(Te-Q,TBDDS2,'LineWidth',2);
% ylabel('Rate per Energy Bin');
% xlabel('Energy - E_{0}^{TT} (eV)');
% legend([a b c],'TT','DT','HT');
% legend boxoff;
% PrettyFigureFormat;
% hold off;


eta_Chi2 = zeros(1,20);
eta_fit  = zeros(1,20);
for i=1:21
    mNu = (i-1)/10;
    Q = RelicNuDebug('Params','KNM1');
    Q.Chi2Twin('Recompute','OFF','Plot','OFF','Syst','OFF','fitPar','E0 Norm Bkg','DeltaChi2',1,'TwinBias_mnuSq',mNu);
    eta_Chi2(i) = Q.etaSensitivity;
    R = MultiRunAnalysis('RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar','eta E0 Norm Bkg',...        % free Parameter!!
        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
        'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
        'fitter','minuit',...                 % minuit standard, matlab to be tried
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...           % final state distribution
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','FSD',...
        'Twin_SameCDFlag','OFF',...
        'Twin_SameIsotopFlag','OFF',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'TwinBias_Q',18573.73,...
        'TwinBias_mnuSq',mNu);

    R.exclDataStart=R.GetexclDataStart(40);
    R.Fit;
    eta_fit(i) = R.FitResult.err(17).*1e10;
end
save('./Results.mat','eta_fit','eta_Chi2');

%load('./RelicNuBkg/UpperLimits/RelicLimit_Twin_BiasmnuSq0_SystOFF_range40_KNM1_mNu E0 Norm Bkg.mat');
%savename = './RelicNuBkg/Chi2Scans/RelicChi2Scan_Twin_BiasmnuSq0_SystOFF_range40_KNM1_[0 5e+11]_mNu E0 Norm Bkg.mat';
%sprintf('%g',eta)
%plotchi2scan(savename);