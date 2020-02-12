%
% Expansion of TBD
%

E0 = 18.574e3; % eV
E  = E0-50:0.1:E0-1;

% Spectrum TBD
S  = @(e,e0,m) (e0-e) .* ((e0-e).^2 - m.^2/2);

% Convolution
 [ zfilt_0 ] = S(E,18574,0);
 [ zfilt_1 ] = gaussfilt(E,S(E,18574,0),0.3);
 [ zfilt_2 ] = gaussfilt(E,S(E,18574,0),0.6);
 [ zfilt_3 ] = gaussfilt(E,S(E,18574,0),0.9);
 zfilt = zfilt_0;
 
% Plot
fig1 = figure(1);
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.7]);
%% Plot neutrino mass
s1=subplot(1,2,1)
% m0=plot(-E0+E,S(E,18574,0)./zfilt ./ S(E(1),18574,0.) .* zfilt(100),'LineWidth',5);
m0=plot(-E0+E,S(E,18574,0)./zfilt ,'LineWidth',5);
hold on
% m1=plot(-E0+E,S(E,18574,0.5)./zfilt_0  ./ S(E(1),18574,0.5) .* zfilt_0(100)  ,'LineWidth',5);
% m2=plot(-E0+E,S(E,18574,0.8)./zfilt_0  ./ S(E(1),18574,0.8) .* zfilt_0(100)  ,'LineWidth',5);
% m3=plot(-E0+E,S(E,18574,1.1)./zfilt_0  ./ S(E(1),18574,1.1) .* zfilt_0(100)  ,'LineWidth',5);
m1=plot(-E0+E,(S(E,18574,0.)./zfilt_1).^-1     ,'LineWidth',5);
m2=plot(-E0+E,(S(E,18574,0.)./zfilt_2).^-1     ,'LineWidth',5);
m3=plot(-E0+E,(S(E,18574,0.)./zfilt_3).^-1    ,'LineWidth',5);
hold off 
grid on
%legend([m0 m1 m2 m3],'m_\nu = 0 eV','m_\nu = 0.2 eV','m_\nu = 0.6 eV','m_\nu = 1 eV','Location','SouthWest');
legend([m0 m1 m2 m3],'0 eV missed-broadening','0.3 eV missed-broadening','0.6 eV missed-broadening','0.9 eV missed-broadening','Location','SouthWest');
ylabel('Spectral Distorsion');
xlabel('E-E_0 (eV)');
PrettyFigureFormat
%ylim([0.995 1]);
%% Plot endpoint
s2=subplot(1,2,2)
%e0=plot(-E0+E,S(E,18574,0)./zfilt ./ S(E(1),18574,0.) .* zfilt(100),'LineWidth',5);
e0=plot(-E0+E,S(E,18574,0)./zfilt  ,'LineWidth',5);
hold on
% e1=plot(-E0+E,S(E,18573.99,0.)./zfilt_0 ./ S(E(1),18573.99,0.) .* zfilt_0(100) ,'LineWidth',5);
% e2=plot(-E0+E,S(E,18573.98,0.)./zfilt_0 ./ S(E(1),18573.98,0.) .* zfilt_0(100) ,'LineWidth',5);
% e3=plot(-E0+E,S(E,18573.97,0.)./zfilt_0 ./ S(E(1),18573.97,0.) .* zfilt_0(100) ,'LineWidth',5);
e1=plot(-E0+E,(S(E,18574,0.)./zfilt_1).^-1   ,'LineWidth',5);
e2=plot(-E0+E,(S(E,18574,0.)./zfilt_2).^-1  ,'LineWidth',5);
e3=plot(-E0+E,(S(E,18574,0.)./zfilt_3).^-1  ,'LineWidth',5);
hold off
grid on
%legend([e0 e1 e2 e3],'\delta E_0 = 0.0 eV','\delta E_0 = -0.01 eV','\delta E_0 = -0.02 eV','\delta E_0 = -0.03 eV','Location','SouthWest');
legend([e0 e1 e2 e3],'0 eV missed-broadening','0.3 eV missed-broadening','0.6 eV missed-broadening','0.9 eV missed-broadening','Location','SouthWest');
ylabel('Spectral Distorsion');
xlabel('E-E_0 (eV)');
%ylim([0.995 1]);
PrettyFigureFormat
linkaxes([s1,s2],'xy');

return;

%% Test Convolution
E0 = 18.574e3; % eV
E  = E0-40:0.1:E0-1;
S  = @(e,e0,m) (e0-e) .* ((e0-e).^2 - m.^2/2);

 figure(999)
 m0=plot(-E0+E,S(E,18574,0)./S(E,18574,0) ,'LineWidth',5);
 hold on
 c1=plot(-E0+E, zfilt_1./S(E,18574,0) ./zfilt_1(300) .* S(E(1),18574,0),'LineWidth',5);
 c2=plot(-E0+E, zfilt_2./S(E,18574,0) ./zfilt_2(300) .* S(E(1),18574,0),'LineWidth',5);
 c3=plot(-E0+E, zfilt_3./S(E,18574,0) ./zfilt_3(300) .* S(E(1),18574,0),'LineWidth',5);
 hold off
 
 