CheckmnuSqDist('Nfit',100);

% Nfit=1000;
% statfluct   = zeros(numel(M.RunData.qU),Nfit);
% statfluct_1=mvnrnd(M.RunData.TBDIS',M.FitCM,Nfit)';
% for j=1:Nfit
%     for i=1:numel(M.RunData.qU)
%         gm=gmdistribution(M.RunData.TBDIS(i),(M.NonPoissonScaleFactor)*M.RunData.TBDIS(i));
%         statfluct(i,j) = random(gm)-M.RunData.TBDIS(i);
%     end
%     statfluct_1(:,j)=statfluct_1(:,j)-M.RunData.TBDIS;
% end

% A=ref_RelicNuBkg_KNM1('TTFSD','OFF','HTFSD','OFF','DTFSD','OFF','mnuSq_i',2);
% B=ref_RelicNuBkg_KNM1('mnuSq_i',2);
% a=semilogy(A.Te-A.Q,A.TBDDS,'LineWidth',2);
% hold on;
% b=semilogy(A.Te-A.Q,A.TBDDS_R,'LineWidth',2);
% c=semilogy(B.Te-B.Q,B.TBDDS,'LineWidth',2);
% d=semilogy(B.Te-B.Q,B.TBDDS_R,'LineWidth',2);
% ylabel('Rate per Energy Bin');
% xlabel('Energy - E_{0} (eV)');
% legend([a b c d],'\beta decay, FSD off','C\nuB, FSD off','\beta decay, FSD on','C\nuB, FSD on');
% legend boxoff;
% xlim([-10 4]);
% ylim([1e-10 1e14]);
% PrettyFigureFormat;

%  load('./RelicNuBkg/Results_mNuFree.mat');
%  %eta_corrected = eta;
%  eta_Chi2_corrected = eta_Chi2;
%  load('./Results_Local_OnlyChi2.mat');
%  eta_local = eta_Chi2;
%  load('./RelicNuBkg/Results_mNuFree_uncorrected.mat');
%  %plot(linspace(0,2,21),eta_fit,'LineWidth',2);
%  hold on;
%  plot(linspace(0,2,21),eta_Chi2,'LineWidth',2);
%  %plot(linspace(0,2,21),eta_corrected,'LineWidth',2);
%  a=plot(linspace(0,2,21),eta_Chi2_corrected,'LineWidth',2);
%  b=plot(linspace(0,2,21),eta_local,'LineWidth',2);
%  ylim([1.39e11 1.95e11]);
%  xlabel('m_{\nu}^{2} (eV^2)');
%  ylabel('1\sigma upper Limit of \eta');
%  %lgd=legend([a b],{'Server','Local'},'Location','northeast','box','off');
%  %lgd.FontSize = 12;
%  PrettyFigureFormat;
%  hold off;

% N=100;
% G=6.67428e-11;
% rho2=linspace(1,11342/2.2,N);  %kg/m³
% rho1=11342; %kg/m³
% R1=3092700; %m
% R2=linspace(1,5400000,N);
% r2=R1+R2+1612e3;
% m1=4*pi./3.*R1^3.*rho1; %kg
% m2=4.*pi./3.*R2.^3.*rho2; %linspace(1,7.117e24,N);
% mu=1-1./(1+sqrt(m2./m1));
% T=2.*pi.*sqrt(r2.^3./(G.*(m1+m2)));
% v_o=(2.*pi.*0.5.*r2)./T;
% for i=1:N
% %     if r2(i)<R1+R2(i)
% %         r2(i)=R1+R2;
% %     end
%     for j=1:N
%         rL(i,j)=r2(j)./(1+sqrt(m2(i)/m1));
%         v(i,j)=sqrt(-2*G*(m1./rL(i,j)+m2(i)./(r2(j)-rL(i,j))-m1./R1-m2(i)./(r2(j)-R1)));
%         if m2(i)>m1.*(r2(j)-R1).^2/R1.^2
%             v(i,j)=0;
%         end
%         vnull(i,j)=sqrt(-2*G*(m1./rL(i,j)-m1./R1));
%     end
% end
% x=linspace(-2*R1,2*R1,300);
% y=linspace(-r2,r2,300);
% vx=[0 0 0];
% for i=1:numel(x)
%     for j=1:numel(y)
%         U(i,j)=-G.*m1./(sqrt(x(i).^2+(y(j)+0.5.*r2).^2))-G.*m2./(sqrt(x(i).^2+(y(j)-0.5.*r2).^2));
%         v_r(i,j)=(2.*pi.*sqrt(x(i)^2+y(j)^2))./T;
%         F(i,j,:)=G.*m1./(x(i).^2+(y(j)+0.5*r2).^2).*[-x(i) -y(j)-r2/2]./sqrt(x(i)^2+(-y(j)-r2/2)^2)+...
%             G.*m2./(x(i).^2+(y(j)-0.5*r2).^2).*[-x(i) -y(j)+r2/2]./sqrt(x(i)^2+(y(j)-r2/2)^2)+...
%             v_r(i,j).^2./sqrt(x(i)^2+y(j)^2).*[-x(i) -y(j)]./sqrt(x(i)^2+y(j)^2)+...
%             2.*[2*pi./T.*vx(2) -vx(1).*2*pi./T];                                  %F(x,y)=[F(x,y,1),F(x,y,2)]
%         if x(i)^2+(y(j)+0.5*r2)^2<R1^2 || (y(j)-0.5*r2)^2+x(i)^2<R2^2
%             U(i,j)=0;
%             F(i,j,:)=0;
%         end
%         %hold on;
%         %plot([x(i) x(i)+1e5*F(i,j,1)],[y(j) y(j)+1e5*F(i,j,2)]);
%     end
% end
% for i=1:10
%     for j=1:10
%         hold on;
%         plot([x(100*i-99) x(100*i-99)+1e5*F(100*i-99,100*j-99,1)],[y(100*j-99) y(100*j-99)+1e5*F(100*i-99,100*j-99,2)]);
%     end
% end
%plot(R2,v);
%hold on;
%plot(R2,vnull);
%plot(R2,r2-R2);
%plot(R2,rL);
%plot([R2(1) R2(end)],[R1 R1]);
%sprintf('Transit velocity: %g m/s',v)
% x=linspace(0,r2,1000);
% y1=sqrt(R1^2-(x-r2).^2);
% y2=sqrt(R2^2-x.^2);
% plot(x,y1);hold on;plot(x,y2);



% D = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
%             'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
%             'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
%             'fixPar','mNu E0 Norm Bkg',...        % free Parameter!!
%             'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
%             'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
%             'minuitOpt','min ; minos',...         % technical fitting options (minuit)
%             'FSDFlag','SibilleFull',...           % final state distribution
%             'ELossFlag','KatrinT2',...            % energy loss function
%             'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
%             'DopplerEffectFlag','FSD',...
%             'Twin_SameCDFlag','OFF',...
%             'Twin_SameIsotopFlag','OFF',...
%             'SynchrotronFlag','ON',...
%             'AngularTFFlag','OFF',...
%             'TwinBias_Q',18573.73,...
%             'TwinBias_mnuSq',0);
% 
%     M = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
%                 'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
%                 'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
%                 'fixPar','E0 Norm Bkg',...    % free Parameter!!
%                 'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
%                 'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
%                 'minuitOpt','min ; minos',...         % technical fitting options (minuit)
%                 'FSDFlag','SibilleFull',...           % final state distribution
%                 'ELossFlag','KatrinT2',...            % energy loss function
%                 'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
%                 'DopplerEffectFlag','FSD',...
%                 'Twin_SameCDFlag','OFF',...
%                 'Twin_SameIsotopFlag','OFF',...
%                 'SynchrotronFlag','ON',...
%                 'AngularTFFlag','OFF',...
%                 'TwinBias_Q',18573.73,...
%                 'TwinBias_mnuSq',-1);
%             
% R = MultiRunAnalysis('RunList','KNM1',...         % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
%                 'chi2','chi2CMShape',...              % uncertainties: statistical or stat + systematic uncertainties
%                 'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
%                 'fixPar','E0 Norm Bkg eta',...    % free Parameter!!
%                 'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
%                 'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
%                 'minuitOpt','min ; minos',...         % technical fitting options (minuit)
%                 'FSDFlag','SibilleFull',...           % final state distribution
%                 'ELossFlag','KatrinT2',...            % energy loss function
%                 'SysBudget',24,...                    % defines syst. uncertainties -> in GetSysErr.m;
%                 'DopplerEffectFlag','FSD',...
%                 'Twin_SameCDFlag','OFF',...
%                 'Twin_SameIsotopFlag','OFF',...
%                 'SynchrotronFlag','ON',...
%                 'AngularTFFlag','OFF',...
%                 'TwinBias_Q',18573.73,...
%                 'TwinBias_mnuSq',-1);
%             
% M.exclDataStart = M.GetexclDataStart(40);
% D.exclDataStart = D.GetexclDataStart(40);
% M.ModelObj.mnuSq_i=4;
% R.ModelObj.mnuSq_i=4;
% R.Fit;
% M.Fit;
% D.Fit;
% 
% %% Global variables
% times = M.ModelObj.qUfrac*M.ModelObj.TimeSec;
% qU    = M.ModelObj.qU; qU    = qU-M.TwinBias_Q; % Energy axis
% NormF = (D.FitResult.par(3+D.ModelObj.nPixels:3+2*D.ModelObj.nPixels-1) + 1)./(M.FitResult.par(3+M.ModelObj.nPixels:3+2*M.ModelObj.nPixels-1) + 1);
% 
% % Spectrum mnu free
% IS  = D.ModelObj.TBDIS; 
% YID = IS./times;
% DIS = D.RunData.TBDIS;
% DIS = DIS./times;
% 
% % Spectrum mnu fixed
% YIM = M.ModelObj.TBDIS; 
% YIM = YIM./times;
% MIS = M.RunData.TBDIS;
% MIS = MIS./times;
% 
% % Spectrum relics
% YIR = R.ModelObj.TBDIS; 
% YIR = YIR./times;
% RIS = R.RunData.TBDIS;
% RIS = RIS./times;
% 
% % Error bar
% err  = (diag(sqrt(D.FitCMShape)));
% err  = err./times;
% err  = err./YIM;
% 
% %% Constraining everything to qULimiteV
% qULimit = -40;
% YID=YID(qU>qULimit);
% DIS=DIS(qU>qULimit);
% MIS=MIS(qU>qULimit);
% YIM=YIM(qU>qULimit);
% RIS=RIS(qU>qULimit);
% YIR=YIR(qU>qULimit);
% err=err(qU>qULimit);
% times=times(qU>qULimit);
% qU=qU(qU>qULimit);
% 
% LocalFontSize = 20;
% 
% fig = figure('Renderer','painters');
% set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);
% 
% plot_title = sprintf('Relic neutrino model: \\eta = %.2g',D.ModelObj.eta);
% prlG = [81 126 102]/255;
% prlB = [50 148 216]/255;
% FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};
% 
% % Ratio
% RSP  = (YIM./DIS);
% RSPd = (YIR./DIS);
% 
% % Plot
% hr1 = plot(qU,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
% hold on;
% hr2 = plot(qU,ones(1,numel(qU)),'color',prlB,'LineWidth',3,'LineStyle',':');
% hr3 = errorbar(qU,(MIS./DIS),err,FitStyleArg{:},'CapSize',0);
% hr4 = plot(qU,RSPd,'color',rgb('Orange'),'LineWidth',3,'LineStyle','-');
% yl2=ylabel('Ratio');
% katrinsim   = sprintf('\\eta=0, m_{\\nu}^2=0 eV^2');
% sterilemod  = sprintf('Fit: \\eta=%.2g, m_{\\nu}^2=4 eV^2 (fix)',R.ModelObj.eta);
% TwinSim     = sprintf('Twin data: m_{\\nu}^2=-1 eV^2');
% norelics    = sprintf('Fit: \\eta=0 (fix), m_{\\nu}^2=4 eV^2 (fix)');
% hl=legend([hr2 hr3 hr1 hr4],{katrinsim,TwinSim,norelics,sterilemod},'Location','southeast','box','off');
% hl.NumColumns=2;
% hl.FontSize = LocalFontSize-2;
% 
% xlim([min(qU-5) 10]);
% %ylim([min((DIS./YI-1)./err) max((DIS./YI-1)./err)]);
% 
% PRLFormat;
% set(gca,'FontSize',LocalFontSize);
% set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
% set(get(gca,'YLabel'),'FontSize',LocalFontSize+4); 
% hl.Position(2) = 0.333;
% hold off;