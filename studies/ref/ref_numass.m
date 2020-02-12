function ul90 = ref_numass(varargin)
%
%            FIT INTEGRAL SPECTRUM
%           Fit Active Neutrino Mass
%    
%                Matlab Fitter
%
%          Th. Lasserre - CEA Saclay
%                February 2017
%

% Initialization
clear par ; clear parnomix;
addpath(genpath('../../Samak2.0'));
set(0,'DefaultFigureWindowStyle','docked')

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('nfit',100,@(x)isfloat(x) && x>=0);
p.addParameter('mnuSq_t',(0.)^2,@(x)isfloat(x));
p.parse(varargin{:});
fign       =    p.Results.fign;
nfit       =    p.Results.nfit;
mnuSq_t    =    p.Results.mnuSq_t;

% Parametrization: True Value
%global A ; 
A=ref_katrin(); A.mnuSq_i=mnuSq_t;
A.ComputeTBDDS(); A.ComputeTBDIS();

% Loop on fits
npar       = 6;
fit_p      = ones(npar,nfit);
fit_xi2min = ones(1,nfit);

tbdisf = @(mb,e0b,bb,nb,qu) A.ComputeTBDISf(qu,mb,e0b,bb,nb,0,0);
tmp = @(mb,e0b,bb,nb,e) interp1(A.qU,tbdisf(mb,e0b,bb,nb,A.qU),e);
myf = fittype(@(mb,e0b,bb,nb,qu) tmp(mb,e0b,bb,nb,qu),...
    'independent','qu','coefficients', {'mb','e0b','bb','nb'});

% Number of Fits
tic
progressbar('...Samak Reference Test - Neutrino Mass Fitting...');
for f=1:1:nfit
    progressbar(f/nfit);
    
    A.ComputeTBDDS(); A.ComputeTBDIS();
    AddStatFluctTBDIS(A);
    Data = [A.qU,A.TBDIS,A.TBDISE];
    
    % Random Init
    i_mnuSq    = (0.2*randn)^2;
    i_E0       = 0.05*randn;
    i_B        = 0.1*A.BKG_RateAllFPDSec*randn;
    i_N        = 1e-3*randn;
    i_m4       = 0;
    i_s4       = 0;
    
    opts = fitoptions('StartPoint',[i_mnuSq,i_E0,i_B,i_N],...
        'Method', 'NonLinearLeastSquares',...
        'Algorithm','Levenberg-Marquardt','Display','off',...
        'Weights',(1./A.TBDIS));
    [fit1,gof,~] = fit(A.qU,A.TBDIS,myf,opts);
    par = [fit1.mb fit1.e0b fit1.bb fit1.nb 0 0];
    ci = confint(fit1,0.68);
    fit_p(:,f)=par; 
    err = [(ci(2,1)-ci(1,1))/2 ...
        (ci(2,2)-ci(1,2))/2 ...
        (ci(2,3)-ci(1,3))/2 ...
        (ci(2,4)-ci(1,4))/2 ...
        0 0];
    chi2min = gof.sse;  % Sum of squares due to error - Equivalent to Chi2
    fit_xi2min(f) = chi2min;
    
    fprintf(2,'-----------------------------------------------------------------------------------\n');
    fprintf('  Processing \t= %g %% \n',f/nfit*100);
    fprintf(2,'-----------------------------------------------------------------------------------\n');
    fprintf(2,'  m^2 \t= %g \t± \t %g \tmeV^2 \t - Init: %g \tmeV^2\n', (A.mnuSq_i+par(1))*1e6,err(1)*1e6,i_mnuSq*1e6);
    fprintf(2,'  dE0 \t= %g \t± \t %g \tmeV \t - Init: %g \tmeV\n', (par(2))*1e3,err(2)*1e3,i_E0*1e3);
    fprintf(2,'  dB \t= %g \t± \t %g \tmcps \t - Init: %g \tmcps\n', (A.BKG_RateSec_i+par(3))*1e3,err(3)*1e3,i_B*1e3);
    fprintf(2,'  dN \t= %g \t± \t %g \t \t - Init: %g \n', par(4),err(4),i_N);
    fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,A.nqU-4);
    fprintf(2,'-----------------------------------------------------------------------------------\n');
    
    figure(10000)
    hdata = errorbar(Data(:,1)-A.Q,...
        Data(:,2)-NuMassModel4par(par,A),...
        Data(:,3),...
        'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
    hold off;
    grid on
    xlabel('qU-E_0 (eV)','FontSize',14);
    ylabel('Residuals','FontSize',14);
    set(gca,'FontSize',12);
    mydata = sprintf('Data: m_{eff}^2=%.3f eV \n',A.mnuSq_i);
    myfit = sprintf('Fit: m^2=%.3f \\pm %.3f eV (\\chi2 / dof=%.1f/%.0f)',sqrt(abs(A.mnuSq_i+par(1))),err(1),chi2min,A.nqU-4);
    legend([hdata hdata],mydata,myfit,'Location','NorthWest') ;
end
toc

%% Plot Results
if nfit<2
figure(fign+5)
fig = figure('Name','Fit and Residuals','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
        
        subplot(2,2,[1 2])
        hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
            'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
        hfit1 = plot(Data(:,1),NuMassModel4par(par,A),...
            'Color','Black','LineWidth',1,'LineStyle','-');
        hfit3 = line([Data(1,1) Data(end,1)],[(A.BKG_RateSec_i+par(3))*A.TimeSec*A.qUfrac(end) (A.BKG_RateSec_i+par(3))*A.TimeSec*A.qUfrac(end)],'LineStyle','--','Color','Red');
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Counts/qU','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        mydata = sprintf('Data: m=%.2f eV - E0=%.2f eV \n',...
            sqrt(mnuSq_t),18575);
        myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n m=%.3f \t\\pm %.3f eV \n E0=%.3f \t\\pm %.3f eV \n B=%.3f \t\\pm %.3f \n N=%.3f \t\\pm %.3f',...
            chi2min,A.nqU-4,sqrt(abs(A.mnuSq_i+par(1))),err(1),(par(2))*1e3,err(2),...
            (A.BKG_RateSec_i+par(3))*1e3,err(3),par(4),err(4));
        mychi2 = sprintf('Data');
        legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
        axis([min(A.qU) max(A.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
        title(sprintf('KATRIN Gaseous Krypton 83m - %s Data and Fit',...
            A.TD),'FontSize',14);
        
        subplot(2,2,[3 4])
        hdata = errorbar(Data(:,1),Data(:,2)-NuMassModel4par(par,A),Data(:,3),...
            'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
        line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Blue');
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Residuals','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        axis([min(A.qU) max(A.qU)+1 min(Data(:,2)-NuMassModel4par(par,A))*2 max(Data(:,2)-NuMassModel4par(par,A))*2]);
        
%         subplot(2,2,4)
%         [h cor] = corplot(err(1:4,1:4));
%         xticklabels({'Position','Width','Amplitude','Offset'})
%         yticklabels({'Position','Width','Amplitude','Offset'})
%         disp(cor);
return;
end

figure(fign)
subplot(2,1,1)
hfit1 = plot(Data(:,1)-A.Q_i,NuMassModel4par(par,A),...
    'LineWidth',1,'LineStyle','-','Color','Black');
hold on;
hfit2 = plot(Data(:,1)-A.Q_i,NuMassModel4par(par,A),...
    'LineWidth',1,'LineStyle','-','Color','Black');
hdata = errorbar(Data(:,1)-A.Q_i,Data(:,2),Data(:,3),...
    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
%errorbar_tick(hdata,200);
hold off
grid on
xlabel('qU-E_0 (eV)','FontSize',14);
ylabel('Counts','FontSize',14);
title(sprintf('KATRIN - %s - %g s - %g mcps',A.TD,A.TimeSec,A.BKG_RateSec_i*1e3));
set(gca,'FontSize',12);
set(gca,'yscale','lin');
mydata = sprintf('Data: m_{eff}=%.2f eV \n',sqrt(abs(A.mnuSq_i)));
myfit = sprintf('Fit: m_{eff}=%.2f \\pm %.2f eV',sqrt(abs(A.mnuSq_i+par(1))),err(1));
mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-4);
legend([hdata hfit1 hfit2],mydata,myfit,mychi2,'Location','NorthEast') ;
axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min(Data(:,2)) max(Data(:,2))*1.1])


subplot(2,1,2)
parnomix = zeros(1,4,1); parnomix = par; parnomix(1) = -A.mnuSq;
hfit = plot(Data(:,1)-A.Q_i,...
    NuMassModel4par(par,A)./NuMassModel4par(parnomix,A)-1,...
    'Color','Black','LineWidth',1,'LineStyle','-');
hold on;
hdata = errorbar(Data(:,1)-A.Q,...
    Data(:,2)./NuMassModel4par(parnomix,A)-1,...
    Data(:,3)./NuMassModel4par(parnomix,A),...
    'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
hold off;
grid on
xlabel('qU-E_0 (eV)','FontSize',14);
ylabel('Spectral distorsion','FontSize',14);
set(gca,'FontSize',12);
a= legend([hdata hfit],'Data/No Mixing - 1','Fit/No Mixing - 1','Location','NorthWest');
axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min((Data(:,2)./NuMassModel4par(parnomix,A)-1)) max(Data(:,3)./NuMassModel4par(parnomix,A))*3])

figure(fign+1)
subplot(2,2,1);
title('Squared Mass','FontSize',12)
nhist(fit_p(1,:)*1e6,'text','pdf','color','sequential');xlabel('meV^2','FontSize',8); grid on
subplot(2,2,2);
title('Endpoint','FontSize',12)
nhist((fit_p(2,:))*1e3,'text','pdf','color','sequential');xlabel('meV','FontSize',8); grid on
subplot(2,2,3);
title('Background','FontSize',12)
nhist(fit_p(3,:)*1e3,'text','pdf','color','sequential');xlabel('mcps','FontSize',8); grid on
subplot(2,2,4);
title('Normalization','FontSize',12)
nhist(((fit_p(4,:))),'text','pdf','color','sequential');xlabel('fraction','FontSize',8); grid on

figure(fign+2)
ndhist(fit_p(1,:)*1e6,((fit_p(2,:)))*1e3);
colorbar
ylabel('Endpoint (meV)','FontSize',14);
xlabel('Mass squared (meV^2)','FontSize',14);
R1 = corrcoef(fit_p(1,:),((fit_p(2,:))));

figure(fign+3)
[t,N,X] = nhist(((fit_p(1,:))*1e6));
% limit on neutrino mass
ul90 = sqrt(X(find(cumsum(N)>0.9.*nfit, 1 )));
fprintf(2,'90 CL Estimate for Active Neutrino Mass Upper Bound: %g meV \n',ul90);
mystd = std(fit_p(1,:)*1e6); %meV^2
fprintf(2,'Variance for Active Neutrino Mass Upper Bound: %g meV \n',sqrt(mystd));
hold on
x = [ul90^2 ul90^2];
y = [0 max(t)];
line(x,y,'Color','red','LineStyle','--')
hold off
ylabel('Bin counts','FontSize',14);
xlabel('Mass squared (meV^2)','FontSize',14);
title(sprintf('KATRIN - %s - %g sec - %g mcps - 90CL = %.2f meV',A.TD,A.TimeSec,A.BKG_RateSec_i*1e3,ul90));
PrettyFigureFormat;

figure(fign+4)
title(sprintf('KATRIN - %s - %g sec - %g mcps',A.TD,A.TimeSec,A.BKG_RateSec_i*1e3));
[t,N,X] = nhist(fit_xi2min);
ylabel('Bin counts','FontSize',14);
xlabel('\chi^2','FontSize',14);
PrettyFigureFormat;

A.DisplayTDBInfo();
end
