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

G=6.67428e-11;
rho2=11342;  %kg/m³
rho1=11342; %kg/m³
R1=3092700; %m
m1=4*pi/3*R1^3*rho1; %kg
r2=R1+R2+1585000;
R2=R1;
m2=4*pi/3*R2.^3*rho2;
mu=1-1./(1+sqrt(m2./m1));
rL=r2./(1+sqrt(m2/m1));
v=sqrt(-2*G*(m1./rL+m2./(r2-rL)-m1./R1-m2./(r2-R1)));
vnull=sqrt(-2*G*(m1./rL-m1./R1));
%plot(R2,v);
%hold on;
%plot(R2,vnull);
%plot(R2,r2-R2);
%plot(R2,rL);
%plot([R2(1) R2(end)],[R1 R1]);
sprintf('Transit velocity: %g m/s',v)
x=linspace(0,r2,1000);
y1=sqrt(R1^2-(x-r2).^2);
y2=sqrt(R2^2-x.^2);
plot(x,y1);hold on;plot(x,y2);