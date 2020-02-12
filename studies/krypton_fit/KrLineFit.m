function [par err chi2min ndof] = KrLineFit(varargin)
%
%            FIT INTEGRAL SPECTRUM
%              Fit a 83mKr Lines
%
%          Th. Lasserre - CEA Saclay
%                January 2018
%

% Initialization
clear par ; clear parnomix;

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('mypub','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Mode','Data',@(x)ismember(x,{'Data','Sim','MoSFitterJuly2017'}));
p.addParameter('CPS','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_Segmentation','PIXEL',@(x)ismember(x,{'OFF','RING','PIXEL'}));
p.addParameter('Pixel',1,@(x)isfloat(x) && x>0);
p.addParameter('nfit',1,@(x)isfloat(x) && x>=0);
p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TD','KrL3_32',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn','Flat20','Flat30','Flat50','Flat60','Flat100','KrK32','KrL3_32','KrL3_32_HS','MoSFitterJuly2017'}));
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','Conv','Voigt','realConv'}));
p.addParameter('SystCorr','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('HVRipples','ON',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('HVRipplesP2PV',0.52,@(x)isfloat(x) && x>0); % V or eV
p.addParameter('TimeYear',10,@(x)isfloat(x) && x>0);

p.parse(varargin{:});
fign              =    p.Results.fign;
mypub             =    p.Results.mypub;
display           =    p.Results.display;
Mode              =    p.Results.Mode;
CPS               =    p.Results.CPS;
FPD_Segmentation  =    p.Results.FPD_Segmentation;
Pixel             =    p.Results.Pixel;
nfit              =    p.Results.nfit;
Fitter            =    p.Results.Fitter;
CovMat            =    p.Results.CovMat;
TD                =    p.Results.TD;
DopplerEffectFlag =    p.Results.DopplerEffectFlag;
SystCorr          =    p.Results.SystCorr;
HVRipples         =    p.Results.HVRipples;
HVRipplesP2PV     =    p.Results.HVRipplesP2PV;
TimeYear          =    p.Results.TimeYear;

clear FitCovErrMatrix; global FitCovErrMatrix;

% Parametrization: True Value
global A ;

switch TD
    case {'KrK32','KrL3_32','KrL3_32_HS'}
        A=InitKrKATRIN_Fit(...
            'TD',TD,'FPD_Pixel',Pixel,...
            'CPS',CPS,'FPD_Segmentation',FPD_Segmentation,...
            'HVRipples',HVRipples,'HVRipplesP2PV',HVRipplesP2PV,...
            'TimeYear',TimeYear,'DopplerEffectFlag',DopplerEffectFlag);
    case 'MoSFitterJuly2017'
        A=InitKrKATRIN_Fit(...
            'TD','MoSFitterJuly2017','FPD_Segmentation','OFF',...
            'HVRipples','OFF','DopplerEffectFlag',DopplerEffectFlag,'CPS','OFF');
end

A.ComputeKrDS(); A.ComputeKrIS();

switch TD
    case 'KrK32'
        LineE_i    = A.K32_E_i;
        LineW_i    = A.K32_W_i;
        LinePhi0_i = A.K32_Phi0_i;
    case {'KrL3_32','KrL3_32_HS'}
        LineE_i    = A.L3_32_E_i;
        LineW_i    = A.L3_32_W_i;
        LinePhi0_i = A.L3_32_Phi0_i;
    case 'MoSFitterJuly2017'
        LineE_i    = A.L3_32_E_i;
        LineW_i    = A.L3_32_W_i;
        LinePhi0_i = A.L3_32_Phi0_i;
        Pixel = 0;
end

switch CovMat
    case 'ON'
        nCMTrials = 10000;
        fName = sprintf('./data/CovMat/WGTSMACE_CovMat_%gTrials_%s_Year_%g.mat',...
            nCMTrials,A.TD,A.TimeYear);
        if exist(fName, 'file') == 2
            CM=importdata(fName);
        else
            %A.ComputeTBDISCovarianceMatrix('nTrials',nCMTrials);
            CM=importdata(fName);
        end
end

% Loop on fits
npar       = 5; fit_p      = ones(npar,nfit); fit_xi2min = ones(1,nfit);

% Init
i_E        = 0; i_W        = 0; i_Bkg      = 0; i_Phi0     = 0;i_N        = 0;

% Data / Sim
switch Mode
    case 'Data'
        [qU krdatamap] = readKrData('TD',TD);
    case 'MoSFitterJuly2017'
        [qU krdatamap] = readKrData('TD',TD);
        A.TimeYear = 1/365.25/86400 ; A.TimeSec = 1; A.qUfrac  = ones(numel(A.qUfrac),1);
        fprintf(2,'MoSFitterJuly2017: Time=%.2f sec\n',A.TimeSec);
        fprintf(2,'MoSFitterJuly2017: Fraction of time in qU bins = %s \n',...
            num2str(A.qUfrac'));
end

% Background / Amplitude Initialization
switch TD
    case 'KrK32'
    case {'KrL3_32'}
        qumin=1;
        tmp = krdatamap(Pixel,1,qumin);
        A.L3_32_Phi0_i  = tmp;
        BKG_i           = krdatamap(Pixel,1,end);
        A.BKG_RateSec_i = BKG_i;
    case {'KrL3_32_HS'}
        qumin=100;
        tmp = krdatamap(Pixel,1,qumin);
        A.L3_32_Phi0_i  = tmp;
        BKG_i           = krdatamap(Pixel,1,end);
        A.BKG_RateSec_i = BKG_i;
end
    
    switch Mode
        case 'Data'
            Count    = (krdatamap(Pixel,1,qumin:end));
            CountErr = (krdatamap(Pixel,2,qumin:end));
            Data = [A.qU(:),Count(:),CountErr(:)];
            A.TimeYear = 1/365.25/86400 ; A.TimeSec = 1;
            A.qUfrac  = ones(numel(A.qUfrac),1);
        case 'Sim'
            A.ComputeKrDS(); A.ComputeKrIS();
            %AddStatFluctKrIS(A);
            %AddStatSystFluctTBDIS(A,'CovMat',CovMat,'SystCorr',SystCorr);
            Data = [A.qU,A.KrIS,A.KrISE];
        case 'MoSFitterJuly2017'
            Count    = (krdatamap(1,1,:));
            CountErr = (krdatamap(1,2,:));
            Data = [A.qU(:),Count(:),CountErr(:)];
            fprintf(2,'MoSFitterJuly2017: Data - Min=%.2f Max=%.2f Ndata=%.0f\n',...
                min(Count(:)), max(Count(:)),numel(Count(:)));
            LinePhi0_i = (Count(1)-Count(end));  A.L3_32_Phi0_i  = LinePhi0_i;
            BKG_i      =  Count(end);            A.BKG_RateSec_i = BKG_i;
            fprintf(2,'MoSFitterJuly2017: L3_32_Phi0_i=%.2f \t BKG_RateSec_i=%.2f \n',...
                A.L3_32_Phi0_i,A.BKG_RateSec_i);
    end
    
    switch display
        case 'ON'
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf(2,'  Initial Parameters \n');
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf(2,'  Init E \t= %.3f eV  \n',   LineE_i);
            fprintf(2,'  Init W \t= %.3f meV \n',   LineW_i);
            fprintf(2,'  Init Phi0 \t= %.3f  \n',   LinePhi0_i);
            fprintf(2,'  Init Bkg \t= %.3f cps \n', BKG_i );
            fprintf(2,'--------------------------------------------------------------\n');
            
    end
  
    
switch Fitter
    case 'Matlab'
        krisf = @(eb,wb,pb,bb,nb,qu) A.ComputeKrISf(qu,eb,wb,pb,bb,nb);
        tmpf = @(eb,wb,pb,bb,nb,e) interp1(A.qU,krisf(eb,wb,pb,bb,nb,A.qU),e);
        myf = fittype(@(eb,wb,pb,bb,nb,qu) tmpf(eb,wb,pb,bb,nb,qu),...
            'independent','qu','coefficients',{'eb','wb','pb','bb'},'problem','nb');
        opts = fitoptions('TolFun',1e-6,'StartPoint',[0 0 0 0],...
            'Method', 'NonLinearLeastSquares',...
            'Display','off',...
            'Weights',(1./A.KrIS));
end

% Number of Fits
tic
progressbar('Krypton Line Fits');
for f=1:1:nfit
    progressbar(f/nfit);
    
    switch Fitter
        case 'Minuit'
            % Initializing Fit
            ParIni = [i_E i_W i_Phi0 i_Bkg i_N];
            parnames = ['E W Phi0 B N'];
            tmparg = sprintf(['set pri -10 ; fix 5 ; set now; min ; imp']);
            Args = {ParIni, Data, '-c',tmparg};
            
            switch CovMat
                case 'OFF'
                    fprintf(2,'No Systematics - Stat. + Pulls\n');
                    [par, err, chi2min, errmat] = fminuit('KrLineChi2',Args{:});
                case 'ON'
                    fprintf(2,'Stat + Systematics (Covariance Matrix)\n');
                    switch SystCorr
                        case 'OFF'
                            FitCovErrMatrix = diag(diag(CM.WGTSMACE_CovMat)) + diag(A.KrIS);
                        case 'ON'
                            FitCovErrMatrix = (CM.WGTSMACE_CovMat) + diag(A.KrIS);
                    end
                    [ par, err, chi2min, errmat] = fminuit('KrLineChi2CovMat',Args{:});
            end
            
            fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
            
        case 'Matlab'
            [fit1,gof,fitinfo] = fit(A.qU,A.KrIS,myf,opts,'problem',0);
            par = [fit1.eb fit1.wb fit1.pb fit1.bb fit1.nb];
            fit_p(:,f)=par; err = [0 0 0 0 0];
            chi2min = gof.sse; fit_xi2min(f) = chi2min;
    end
    
    % Definition of Amplitude / Background
    E_Fit             = (LineE_i+par(1));
    E_FitError        = (err(1));
    
    W_Fit             = (LineW_i+par(2))*1e3;
    W_FitError        = (par(2))*1e3;
    
    Phi0_Fit          = (A.L3_32_Phi0_i+par(3))*(A.TimeSec.*A.qUfrac(1));
    Phi0_FitError     = err(3)/(A.TimeSec.*A.qUfrac(1));
    
    Offset_Fit        = (A.BKG_RateSec_i+par(4))*(A.TimeSec.*A.qUfrac(1));
    Offset_FitError   = err(4)/(A.TimeSec.*A.qUfrac(1));
    
    switch display
        case 'ON'
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf('  Processing \t= %g %% \n',f/nfit*100);
            fprintf(2,'  E \t= %.3f \t± \t %g \t eV \n',E_Fit,E_FitError);
            fprintf(2,'  W \t= %.3f \t± \t %g \t meV \n',W_Fit,W_FitError);
            fprintf(2,'  Phi0 \t= %.3f \t± \t %g \t \n',Phi0_Fit,Phi0_FitError);
            fprintf(2,'  Bkg \t= %.3f \t± \t %g \t cps \n',Offset_Fit,Offset_FitError);
            fprintf(2,'  N \t= %g \t± \t %g \t fixed \n',(par(5)),err(5));
            ndof = A.nqU-4;fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
            fprintf(2,'----------------------------------------------------\n');
    end
end
toc

%% Plot Results
str = sprintf('Fit and Residuals - Pixel %.0f',Pixel);
fig = figure('Name',str,'NumberTitle','off','rend','painters','pos',[10 10 1400 800]);

subplot(2,2,1)
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
    'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
hfit1 = plot(Data(:,1),KrLineModel4par(par),...
    'Color','Black','LineWidth',1,'LineStyle','-');
hfit3 = line([Data(1,1) Data(end,1)],[par(4) par(4) ],'LineStyle','--','Color','Red');
hold off
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Counts/qU','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
    LineE_i,LineW_i);
myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n E=%.3f \t\\pm %.3f eV \n W=%.0f \t\\pm %.0f meV \n A=%.3f \t\\pm %.3f \n O=%.3f \t\\pm %.3f',...
    chi2min,A.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError,...
    Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
mychi2 = sprintf('Data');
legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
axis([min(A.qU) max(A.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
title(sprintf('KATRIN Gaseous Krypton 83m - {%s} Data and Fit',...
    A.TD),'FontSize',14);

subplot(2,2,3)
hdata = errorbar(Data(:,1),Data(:,2)-KrLineModel4par(par),Data(:,3),...
    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Blue');
hold off
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Residuals','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
axis([min(A.qU) max(A.qU)+1 min(Data(:,2)-KrLineModel4par(par))*2 max(Data(:,2)-KrLineModel4par(par))*2]);

subplot(2,2,2)
hfit2 = plot(A.Te,(KrLineModelDiff4par(par)),...
    'Color','Black','LineWidth',1,'LineStyle','-');
grid on
xlabel('Te (eV)','FontSize',14);
ylabel('Counts/Te','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
title(sprintf('Differential Spectrum - Doppler %s (fixed \\sigma=%.1f eV)',...
    A.DopplerEffectFlag,A.DE_sigma),'FontSize',14);

subplot(2,2,4)
%figure('Name','Fit Correlation Matrix','NumberTitle','off','rend','painters','pos',[10 10 500 500])
%title(sprintf('KATRIN Krypton 83m Fit - Correlation Matrix - %s - Doppler %s\n',A.TD,A.DopplerEffectFlag));
[h cor] = corrplot(errmat(1:4,1:4));
xticklabels({'Position','Width','Amplitude','Offset'})
yticklabels({'Position','Width','Amplitude','Offset'})
disp(cor);

        
        %             figure(fign+1000)
        %             [h cor] = corplot(errmat(1:4,1:4));
        %             xticklabels({'Position','Width','Amplitude','Background'})
        %             yticklabels({'Position','Width','Amplitude','Background'})
        %             title(sprintf('KATRIN Krypton - Correlation Matrix - %s - Doppler %s\n',A.TD,A.DopplerEffectFlag));
        %             fprintf(2,'Error Matrixï¿½\n');
        %             disp(cor);
        
        switch Mode
            case 'MoSFitterJuly2017'
                switch pub
                    case 'YES'
                        myname = sprintf('./figures/MoSFitterJuly2017_2.eps');
                        fprintf(2,'publish: %s\n',myname); %publish_figure(1000+fign,myname);
                end
            case {'DATA','Sim'}
                myname = sprintf('./figures/krypton_%s_pixel%d_1.eps',...
                    TD,Pixel);
                fprintf(2,'publish: %s\n',myname);publish_figure(fig,myname);
        end
end


if nfit<5
    % For Saving:
    par(1) = LineE_i         + par(1);
    par(2) = LineW_i         + par(2);
    par(3) = A.L3_32_Phi0_i  + par(3);
    par(4) = A.BKG_RateSec_i + par(4);
end

 switch Mode
    case 'Data'
       return;
    case 'MoSFitterJuly2017'
       return;
 end


    figure(fign)
    subplot(2,1,1)
    hfit1 = plot(Data(:,1),KrLineModel4par(par),...
        'LineWidth',1,'LineStyle','-','Color','Black');
    hold on;
    hfit2 = plot(Data(:,1),KrLineModel4par(par),...
        'LineWidth',1,'LineStyle','-','Color','Black');
    hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
        'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
    %errorbar_tick(hdata,200);
    hold off
    grid on
    xlabel('qU (eV)','FontSize',10);
    ylabel('Counts','FontSize',10);
    title(sprintf('KATRIN Integral - %s - %g y - %g mcps',A.TD,A.TimeYear,A.BKG_RateSec_i));
    set(gca,'FontSize',12);
    set(gca,'yscale','lin');
    mydata = sprintf('Data: E=%.2f eV \n',LineE_i);
    myfit = sprintf('Fit: E=%.2f \\pm %.2f eV',LineE_i+par(1),err(1));
    mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-4);
    legend([hdata hfit1 hfit2],mydata,myfit,mychi2,'Location','NorthEast') ;% legend(a,'boxoff');
    axis([min(A.qU) max(A.qU)+1 min(Data(:,2)) max(Data(:,2))*1.1])
    
    subplot(2,1,2)
    parnomix = zeros(1,4,1); parnomix = par; parnomix(1) = 0; parnomix(2) = 0;
    hfit = plot(Data(:,1),...
        KrLineModel4par(par)./KrLineModel4par(parnomix)-1,...
        'Color','Black','LineWidth',1,'LineStyle','-');
    hold on;
    hdata = errorbar(Data(:,1),...
        Data(:,2)./KrLineModel4par(parnomix)-1,...
        Data(:,3)./KrLineModel4par(parnomix),...
        'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
    %errorbar_tick(hdata,200);
    hold off;
    grid on
    xlabel('qU (eV)','FontSize',10);
    ylabel('Spectral distorsion','FontSize',10);
    set(gca,'FontSize',12);
    a= legend([hdata hfit],'Data/No Mixing - 1','Fit/No Mixing - 1','Location','NorthWest');
    %legend(a,'boxoff');
    %axis([min(A.qU) max(A.qU)+1 min((Data(:,2)./KrLineModel4par(parnomix)-1)) max(Data(:,3)./KrLineModel4par(parnomix))*3])

    switch mypub
        case 'ON'
           myname = sprintf('./figures/f1-krkatrin_%s_pixel%.0f.eps',...
                TD,Pixel); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign,myname);
    end
    
    figure(fign+1)
    subplot(2,2,1);
    title('K32_E','FontSize',10)
    nhist(fit_p(1,:)*1e3,'text','pdf','color','sequential');xlabel('meV','FontSize',8); grid on
    subplot(2,2,2);
    title('K32_W','FontSize',12)
    nhist((fit_p(2,:))*1e3,'text','pdf','color','sequential');xlabel('meV','FontSize',8); grid on
    subplot(2,2,3);
    title('Normalzation','FontSize',12)
    nhist(fit_p(3,:),'text','pdf','color','sequential');xlabel('no unit','FontSize',8); grid on
    subplot(2,2,4);
    title('Background','FontSize',12)
    nhist(((fit_p(4,:))),'text','pdf','color','sequential');xlabel('mcps','FontSize',8); grid on
    switch mypub
        case 'ON'
           myname = sprintf('./figures/f2-krkatrin_%g-mcps_%g-numsq.eps',...
                A.BKG_RateAllFPDSec*1e3,0); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign+1,myname);
    end
    
    figure(fign+2)
    ndhist(fit_p(1,:)*1e3,((fit_p(2,:)))*1e3);
    colorbar
    ylabel('E (meV)','FontSize',10);
    xlabel('W (meV)','FontSize',10);
    R1 = corrcoef(fit_p(1,:),((fit_p(2,:))));
    switch mypub
        case 'ON'
           myname = sprintf('./figures/f3-krkatrin_%g-mcps_%g-numsq.eps',...
                A.BKG_RateAllFPDSec*1e3,0); 
            fprintf(2,'publish: %s\n',myname);publish_figure(fign+2,myname);
    end
         
end
