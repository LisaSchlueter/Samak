function fit_p = KrLineFitSim(varargin)
%
%            FIT INTEGRAL SPECTRUM
%          Fit a 83mKr Simulated Lines
%
%
%          Th. Lasserre - CEA Saclay
%                January 2018
%

% Initialization
clear par ; clear parnomix;
format long

% Parser
p = inputParser;
p.addParameter('fign',1000,@(x)isfloat(x) && x>0);
p.addParameter('mypub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Pixel',1,@(x)isfloat(x) && x>0);
p.addParameter('nfit',1,@(x)isfloat(x) && x>=0);
p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
p.addParameter('TD','DummyKrdc1',@(x)ismember(x,{'DummyKrdc1','KrK32','KrL3_32','KrL3_32_HS'}));
p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('SystCorr','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('DCKr1','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});
mypub             =    p.Results.mypub;
display           =    p.Results.display;
Pixel             =    p.Results.Pixel;
nfit              =    p.Results.nfit;
Fitter            =    p.Results.Fitter;
TD                =    p.Results.TD;
CovMat            =    p.Results.CovMat;
SystCorr          =    p.Results.SystCorr;
DCKr1             =    p.Results.DCKr1;

clear FitCovErrMatrix; global FitCovErrMatrix;

% Parametrization: True Value
global B ;  %Object (Modell) for Simulation ---->MACE
global A;   %Object (Modell) for Fit used in KrLineModel...m ---->MACER
switch TD
    case {'KrK32','KrL3_32','KrL3_32_HS','DummyKrdc1'}
        B=InitKrKATRIN_DummyKrdc1(...
            'TD',TD,'FPD_Pixel',Pixel);
end

B.ComputeKrDS(); B.ComputeKrIS();
BKG_i          = B.BKG_RateSec_i;

switch TD
    case 'KrK32'
        LineE_i    = B.K32_E_i;
        LineW_i    = B.K32_W_i;
        LinePhi0_i = B.K32_Phi0_i;
    case {'KrL3_32','KrL3_32_HS','DummyKrdc1'}
        LineE_i    = B.L3_32_E_i;
        LineW_i    = B.L3_32_W_i;
        LinePhi0_i = B.L3_32_Phi0_i;
end

switch CovMat
    case 'ON'
        nCMTrials = 10000;
        fName = sprintf('./data/CovMat/WGTSMACE_CovMat_%gTrials_%s_Year_%g.mat',...
            nCMTrials,B.TD,B.TimeYear);
        if exist(fName, 'file') == 2
            CM=importdata(fName);
        else
            %A.ComputeTBDISCovarianceMatrix('nTrials',nCMTrials);
            CM=importdata(fName);
        end
end

% Loop on fits
npar       = 5; fit_p      = ones(npar,nfit); fit_xi2min = ones(1,nfit);
ndof       = B.nqU-4; %1 par is fixed

% Init
i_E        = 0; i_W        = 0; i_Bkg      = 0; i_Phi0     = 0;i_N        = 0;

% Generate Data
B.ComputeKrDS(); B.ComputeKrIS();

switch display
    case 'ON'
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  Initial Parameters \n');
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  Init E \t= %.3f \t eV  \n',      LineE_i);
        fprintf(2,'  Init W \t= %.3f \t meV \n',      LineW_i);
        switch B.CPS
            case 'OFF'
                fprintf(2,'  Init Phi0 \t= %.3f \t counts\n', LinePhi0_i);
                fprintf(2,'  Init Bkg \t= %.3f \t counts \n', BKG_i );
            case 'ON'
                fprintf(2,'  Init Phi0 \t= %.3f \t cps\n',    LinePhi0_i);
                fprintf(2,'  Init Bkg \t= %.3f \t cps \n',    BKG_i );
        end
        fprintf(2,'--------------------------------------------------------------\n');
        
end

switch Fitter
    case 'Matlab'
        krisf = @(eb,wb,pb,bb,nb,qu) B.ComputeKrISf(qu,eb,wb,pb,bb,nb);
        tmpf = @(eb,wb,pb,bb,nb,e) interp1(B.qU,krisf(eb,wb,pb,bb,nb,B.qU),e);
        myf = fittype(@(eb,wb,pb,bb,nb,qu) tmpf(eb,wb,pb,bb,nb,qu),...
            'independent','qu','coefficients',{'eb','wb','pb','bb'},'problem','nb');
        opts = fitoptions('TolFun',1e-6,'StartPoint',[0 0 0 0],...
            'Method', 'NonLinearLeastSquares',...
            'Display','off',...
            'Weights',(1./B.KrIS));
end

% Number of Fits
tic
progressbar('Krypton Line Fits');
for f=1:1:nfit
    progressbar(f/nfit);
    switch DCKr1
        case 'OFF'
            A.ComputeKrDS(); A.ComputeKrIS();
            AddStatFluctKrIS(A);
            %AddStatSystFluctTBDIS(A,'CovMat',CovMat,'SystCorr',SystCorr);
            Data = [A.qU,A.KrIS,A.KrISE];
        case 'ON'
            nfit = 1;
            file=importdata('kdc_v1_samak.txt');
            Data = [A.qU,file(:,2),file(:,3)];
    end
   
    switch Fitter
        case 'Minuit'
            % Initializing Fit
            ParIni = [i_E i_W i_Phi0 i_Bkg i_N];
            parnames = ['E W Phi0 B N'];
            tmparg = sprintf(['set pri -10 ; fix 5 ; set now; min ; imp; set pri -10']);
            Args = {ParIni, Data, '-c',tmparg};
            
            switch CovMat
                case 'OFF'
                    fprintf(2,'No Systematics - Stat. + Pulls\n');
                    [par, err, chi2min, errmat] = fminuit('KrLineChi2',Args{:});
                case 'ON'
                    fprintf(2,'Stat + Systematics (Covariance Matrix)\n');
                    switch SystCorr
                        case 'OFF'
                            FitCovErrMatrix = diag(diag(CM.WGTSMACE_CovMat)) + diag(B.KrIS);
                        case 'ON'
                            FitCovErrMatrix = (CM.WGTSMACE_CovMat) + diag(B.KrIS);
                    end
                    [ par, err, chi2min, errmat] = fminuit('KrLineChi2CovMat',Args{:});
            end
            
            fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
            
        case 'Matlab'
            [fit1,gof,fitinfo] = fit(B.qU,B.KrIS,myf,opts,'problem',0);
            par = [fit1.eb fit1.wb fit1.pb fit1.bb fit1.nb];
            fit_p(:,f)=par; err = [0 0 0 0 0];
            chi2min = gof.sse; fit_xi2min(f) = chi2min;
    end
    
    % Definition of Amplitude / Background    
    E_Fit             = (LineE_i+par(1));
    E_FitError        = (err(1));
    
    W_Fit             = (LineW_i+par(2))*1e3;
    W_FitError        = (err(2))*1e3;
 
    switch B.CPS
        case 'OFF'
            Phi0_Fit          = (B.L3_32_Phi0_i+par(3))*(B.TimeSec.*B.qUfrac(1));
            Phi0_FitError     = err(3)*(B.TimeSec.*B.qUfrac(1));
            Offset_Fit        = (B.BKG_RateSec_i+par(4))*(B.TimeSec.*B.qUfrac(1));
            Offset_FitError   = err(4)*(B.TimeSec.*B.qUfrac(1));
        case 'ON'
            Phi0_Fit          = (B.L3_32_Phi0_i+par(3));
            Phi0_FitError     = err(3);
            Offset_Fit        = (B.BKG_RateSec_i+par(4));
            Offset_FitError   = err(4);
    end
            
    switch display
        case 'ON'
            
            B.DisplayKrInfo;
            fprintf(2,'----------------------------------------------------\n');
            fprintf(2,'Fit Kr83m Line \n');
            fprintf(2,'----------------------------------------------------\n');
            fprintf(2,'  E \t= %.3f \t� \t %g \t eV \n',E_Fit,E_FitError);
            fprintf(2,'  W \t= %.3f \t� \t %g \t meV \n',W_Fit,W_FitError);
            fprintf(2,'  Phi0 \t= %.3f \t� \t %g \t \n',Phi0_Fit,Phi0_FitError);
            fprintf(2,'  Bkg \t= %.3f \t� \t %g \t  \n',Offset_Fit,Offset_FitError);
            fprintf(2,'  N \t= %g \t� \t %g \t fixed \n',(par(5)),err(5));
            fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
            fprintf(2,'----------------------------------------------------\n');
    end
end
toc

%% Plot Results
str1 = sprintf('Fit and Residuals - Pixel %.0f',Pixel);
fig1 = figure('Name',str1,'NumberTitle','off','rend','painters','pos',[10 10 1400 800]*0.9);
subplot(2,1,1)
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
    'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
hfit1 = plot(Data(:,1),KrLineModel4par(par),...
    'Color','Black','LineWidth',1,'LineStyle','-');
hfit3 = line([Data(1,1) Data(end,1)],[Offset_Fit Offset_Fit],'LineStyle','--','Color','Red');
hold off
grid on
xlabel('qU (eV)','FontSize',14);
switch B.CPS
    case 'OFF'
        ylabel('Counts','FontSize',14);
    case 'ON'
        ylabel('Rates','FontSize',14);
end
    
set(gca,'FontSize',14);
set(gca,'yscale','lin');
mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
    LineE_i,LineW_i);
myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n E=%.3f \t\\pm %.3f eV \n W=%.0f \t\\pm %.0f meV \n A=%.3f \t\\pm %.3f \n O=%.3f \t\\pm %.3f',...
    chi2min,B.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError,...
    Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
mychi2 = sprintf('Data');
legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
axis([min(B.qU) max(B.qU) 0.*min(Data(:,2)) max(Data(:,2))*1.2])
title(sprintf('Systematic Effect: TF - {%s} Simulation (non-rel. TF) and Fit (Model: rel TF)',...
    B.TD),'FontSize',14);

subplot(2,1,2)
hdata = errorbar(Data(:,1),Data(:,2)-KrLineModel4par(par),Data(:,3),...
    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
line([min(B.qU) max(B.qU)],[0 0 ],'LineStyle','--','Color','Blue');
hold off
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Residuals','FontSize',14);
axis([min(B.qU) max(B.qU) -0.1 0.1]);
set(gca,'FontSize',14);
set(gca,'yscale','lin');

subplot(2,2,2)
hfit2 = plot(A.Te,(KrLineModelDiff4par(par)),...
    'Color','Black','LineWidth',1,'LineStyle','-');
grid on
xlabel('Te (eV)','FontSize',14);
ylabel('Counts/Te','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
title(sprintf('Differential Spectrum - Doppler %s (fixed \\sigma=%.3f eV)',...
    A.DopplerEffectFlag,A.DE_sigma),'FontSize',14);

% subplot(2,2,4)
% %figure('Name','Fit Correlation Matrix','NumberTitle','off','rend','painters','pos',[10 10 500 500])
% %title(sprintf('KATRIN Krypton 83m Fit - Correlation Matrix - %s - Doppler %s\n',A.TD,A.DopplerEffectFlag));
% [h cor] = corrplot(errmat(1:4,1:4));
% xticklabels({'Position','Width','Amplitude','Offset'})
% yticklabels({'Position','Width','Amplitude','Offset'})
% disp(cor);

switch mypub
    case 'ON'
        myname = sprintf('./figures/krypton_%s_pixel%d_1.eps',...
            TD,Pixel);
        fprintf(2,'publish: %s\n',myname);publish_figure(fig1,myname);
end


if nfit<2
    % For Saving:
    par(1) = LineE_i         + par(1);
    par(2) = LineW_i         + par(2);
    par(3) = B.L3_32_Phi0_i  + par(3);
    par(4) = B.BKG_RateSec_i + par(4);
return;
end

str2 = sprintf('Distribution of Fit Parameters - %.0f fits',nfit);
fig2 = figure('Name',str2,'NumberTitle','off','rend','painters','pos',[30 30 1400 800]*0.9);
hbw = 20;
subplot(2,2,1);
%nhist(fit_p(1,:)*1e3,'text','pdf','color','summer','binfactor',1);
h=histfit(fit_p(1,:)*1e3,hbw,'normal');h(1).FaceColor = [.7 .7 1];
pd = fitdist(fit_p(1,:)'*1e3,'Normal');
pdpar = sprintf('mean = %.1f \n \\sigma = %.1f',pd.mu,pd.sigma);
legend(h,pdpar,'Location','NorthEast') ;
xlabel('meV','FontSize',12); ylabel('Counts','FontSize',12); grid on
title('Line Position','FontSize',12)
subplot(2,2,2);
%nhist((fit_p(2,:))*1e3,'text','pdf','color','summer','binfactor',1);
h=histfit(fit_p(2,:)*1e3,hbw,'normal');h(1).FaceColor = [.7 .7 1];
pd = fitdist(fit_p(2,:)'*1e3,'Normal');
pdpar = sprintf('mean = %.1f \n \\sigma = %.1f',pd.mu,pd.sigma);
legend(h,pdpar,'Location','NorthEast') ;
xlabel('meV','FontSize',12); ylabel('Counts','FontSize',12); grid on
title('Line Width','FontSize',12)
subplot(2,2,3);
%nhist(fit_p(3,:),'text','pdf','color','summer','binfactor',1);
h=histfit(fit_p(3,:),hbw,'normal'); ylabel('Counts','FontSize',12); h(1).FaceColor = [.7 .7 1];
pd = fitdist(fit_p(3,:)','Normal');
pdpar = sprintf('mean = %.1f \n \\sigma = %.1f',pd.mu,pd.sigma);
legend(h,pdpar,'Location','NorthEast') ;
xlabel('cps','FontSize',12); grid on
title('Line Amplitude','FontSize',12)
subplot(2,2,4);
%nhist(((fit_p(4,:))),'text','pdf','color','summer','binfactor',1);
h=histfit(fit_p(4,:),hbw,'normal'); ylabel('Counts','FontSize',12) ; h(1).FaceColor = [.7 .7 1];
pd = fitdist(fit_p(4,:)','Normal');
pdpar = sprintf('mean = %.1f \n \\sigma = %.1f',pd.mu,pd.sigma);
legend(h,pdpar,'Location','NorthEast') ;
xlabel('cps','FontSize',12); grid on
title('Offset','FontSize',12)
switch mypub
    case 'ON'
        myname = sprintf('./figures/krypton_%s_pixel%d_2.eps',...
            TD,Pixel);
        fprintf(2,'publish: %s\n',myname);publish_figure(fig2,myname);
end

str3 = sprintf('\\chi2 - %.0f fits',nfit);
fig3 = figure('Name',str3,'NumberTitle','off','rend','painters','pos',[300 300 600 600]*0.9);
title('\chi^2 Distribution','FontSize',12)
[t, N, X] = nhist(fit_xi2min,'text','pdf','color','color',{[.7 .7 1]},'binfactor',1);
hold on
plot(X,smooth(chi2pdf(X,ndof)),'LineWidth',2,'LineStyle','--');
hold off
xlabel('\chi2','FontSize',12); grid on
switch mypub
    case 'ON'
        myname = sprintf('./figures/krypton_%s_pixel%d_3.eps',...
            TD,Pixel);
        fprintf(2,'publish: %s\n',myname);publish_figure(fig3,myname);
end

cm=fit_p(1:4,:)';
corrplotm(cm,'rows','complete','varNames',{'E','W','A','O'},'testR','on');
switch mypub
    case 'ON'
        myname = sprintf('./figures/krypton_%s_pixel%d_4.eps',...
            TD,Pixel);
        fprintf(2,'publish: %s\n',myname);publish_figure(4,myname);
end
fit_p=cm;
end
