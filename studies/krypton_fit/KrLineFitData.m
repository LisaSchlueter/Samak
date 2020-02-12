function [par, err, fit_xi2min, ndof] = KrLineFitData(varargin)
%
% FIT INTEGRAL SPECTRUM DATA
% 83mKr Lines K32 and L3-32
% Written for the KATRIN Krypton First Lights
%
% KrLineFit('TD','KrL3_32',...
%         'mypub','OFF',...
%         'FPD_Segmentation','PIXEL','Pixel',1,...
%         'DopplerEffectFlag','ON',...
%         'HVRipples','ON','HVRipplesP2PV',0.52);
%
% Th. Lasserre - CEA Saclay
% Last Updated: January 2018
%

% Initialization
clear par ; clear parnomix;

% Parser
p = inputParser;
p.addParameter('mypub','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('CPS','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_Segmentation','PIXEL',@(x)ismember(x,{'OFF','RING','PIXEL'}));
p.addParameter('Pixel',1,@(x)isfloat(x) && x>0);
p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TD','DummyKrdc1',@(x)ismember(x,{'KrK32','KrL3_32','KrL3_32_HS','DummyKrdc1'}));
p.addParameter('DopplerEffectFlag','Conv',@(x)ismember(x,{'OFF','Conv','Voigt','realConv'}));
p.addParameter('SystCorr','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('HVRipples','ON',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('HVRipplesP2PV',0.44,@(x)isfloat(x) && x>0); % V or eV
p.addParameter('KrObject','',@(x) isa(x, 'Kr')); 
p.parse(varargin{:});

mypub             =    p.Results.mypub;
display           =    p.Results.display;
CPS               =    p.Results.CPS;
FPD_Segmentation  =    p.Results.FPD_Segmentation;
Pixel             =    p.Results.Pixel;
Fitter            =    p.Results.Fitter;
CovMat            =    p.Results.CovMat;
TD                =    p.Results.TD;
DopplerEffectFlag =    p.Results.DopplerEffectFlag;
SystCorr          =    p.Results.SystCorr;
HVRipples         =    p.Results.HVRipples;
HVRipplesP2PV     =    p.Results.HVRipplesP2PV;
A                 =    p.Results.KrObject;

clear FitCovErrMatrix; global FitCovErrMatrix;

% Parametrization: True Value
%global A ;
switch TD
    case {'KrK32','KrL3_32','KrL3_32_HS'}
        A=InitKrKATRIN_Fit(...
            'TD',TD,'FPD_Pixel',Pixel,...
            'DopplerEffectFlag',DopplerEffectFlag,...
            'CPS',CPS,...
            'FPD_Segmentation',FPD_Segmentation,...
            'HVRipples',HVRipples,'HVRipplesP2PV',HVRipplesP2PV);
         case {'DummyKrdc1'}
        A=InitKrKATRIN_DummyKrdc1();  % Benchmark setting for Kr data challenge (2018)
end
A.ComputeKrDS(); A.ComputeKrIS();

switch TD
    case 'KrK32'
        LineE_i    = A.K32_E_i;
        LineW_i    = A.K32_W_i;
        LinePhi0_i = A.K32_Phi0_i;
    case {'KrL3_32','KrL3_32_HS','DummyKrdc1'}
        LineE_i    = A.L3_32_E_i;
        LineW_i    = A.L3_32_W_i;
        LinePhi0_i = A.L3_32_Phi0_i;
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
npar       = 5; fit_p      = ones(npar,1); fit_xi2min = ones(1,1);

% Init
i_E        = 0; i_W        = 0; i_Bkg      = 0; i_Phi0     = 0;i_N        = 0;

% Read Data File
[~, krdatamap] = readKrData('TD',TD);

% Background / Amplitude Initialization / Fill Data
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
   case {'DummyKrdc1'}
       qumin=1;
        tmp = krdatamap(1,1,qumin);
        A.L3_32_Phi0_i  = tmp;
        BKG_i           = krdatamap(1,1,end);
        A.BKG_RateSec_i = BKG_i; 
end
Count    = (krdatamap(Pixel,1,qumin:end));
CountErr = (krdatamap(Pixel,2,qumin:end));
Data = [A.qU(:),Count(:),CountErr(:)];
A.TimeYear = 1/365.25/86400 ; A.TimeSec = 1;
A.qUfrac  = ones(numel(A.qUfrac),1);

switch display
    case 'ON'
        
        A.DisplayKrInfo;

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

switch Fitter
    case 'Minuit'
        % Initializing Fit
        ParIni = [i_E i_W i_Phi0 i_Bkg i_N];
        parnames = ['E W Phi0 B N'];
        tmparg = sprintf(['set pri -10 ; fix 5 ; set now; min ; imp']);
        DataKr = {Data, A};
        Args = {ParIni, DataKr, '-c',tmparg};
        
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
        
        fit_p(:)=par;  fit_xi2min = chi2min;
        
    case 'Matlab'
        [fit1,gof,~] = fit(A.qU,A.KrIS,myf,opts,'problem',0);
        par = [fit1.eb fit1.wb fit1.pb fit1.bb fit1.nb];
        fit_p(:)=par; err = [0 0 0 0 0];
        chi2min = gof.sse; fit_xi2min(f) = chi2min;
end

% Definition of Amplitude / Background
E_Fit             = (LineE_i+par(1));
E_FitError        = (err(1));

W_Fit             = (LineW_i+par(2))*1e3;
W_FitError        = (err(2))*1e3;

Phi0_Fit          = (A.L3_32_Phi0_i+par(3))*(A.TimeSec.*A.qUfrac(1));
Phi0_FitError     = err(3)/(A.TimeSec.*A.qUfrac(1));

Offset_Fit        = (A.BKG_RateSec_i+par(4))*(A.TimeSec.*A.qUfrac(1));
Offset_FitError   = err(4)/(A.TimeSec.*A.qUfrac(1));

switch display
    case 'ON'
        fprintf(2,'----------------------------------------------------\n');
        fprintf(2,'Fit Kr83m L3-32 Line \n');
        fprintf(2,'----------------------------------------------------\n');
        fprintf(2,'  E \t= %.3f +- \t %g \t eV \n',E_Fit,E_FitError);
        fprintf(2,'  W \t= %.3f +- \t %g \t meV \n',W_Fit,W_FitError);
        fprintf(2,'  Phi0 \t= %.3f +- \t %g \t \n',Phi0_Fit,Phi0_FitError);
        fprintf(2,'  Bkg \t= %.3f +- \t %g \t cps \n',Offset_Fit,Offset_FitError);
        fprintf(2,'  N \t= %g +- \t %g \t fixed \n',(par(5)),err(5));
        ndof = A.nqU-4;fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
        fprintf(2,'----------------------------------------------------\n');
end

% %% Plot Results
% str = sprintf('Fit and Residuals - Pixel %.0f',Pixel);
% fig = figure('Name',str,'NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
% 
% subplot(2,2,1)
% hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
%     'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
% hold on
% hfit1 = plot(Data(:,1),KrLineModel4par(par),...
%     'Color','Black','LineWidth',1,'LineStyle','-');
% hfit3 = line([Data(1,1) Data(end,1)],[Offset_Fit Offset_Fit],'LineStyle','--','Color','Red');
% hold off
% grid on
% xlabel('qU (eV)','FontSize',14);
% ylabel('Counts/qU','FontSize',14);
% set(gca,'FontSize',14);
% set(gca,'yscale','lin');
% mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
%     LineE_i,LineW_i);
% myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n E=%.3f \t\\pm %.3f eV \n W=%.0f \t\\pm %.0f meV \n A=%.3f \t\\pm %.3f \n O=%.3f \t\\pm %.3f',...
%     chi2min,A.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError,...
%     Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
% mychi2 = sprintf('Data');
% legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
% axis([min(A.qU) max(A.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
% title(sprintf('KATRIN Gaseous Krypton 83m - {%s} Data and Fit',...
%     A.TD),'FontSize',14);
% 
% subplot(2,2,3)
% hdata = errorbar(Data(:,1),Data(:,2)-KrLineModel4par(par),Data(:,3),...
%     'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
% hold on
% line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Blue');
% hold off
% grid on
% xlabel('qU (eV)','FontSize',14);
% ylabel('Residuals','FontSize',14);
% set(gca,'FontSize',14);
% set(gca,'yscale','lin');
% axis([min(A.qU) max(A.qU)+1 min(Data(:,2)-KrLineModel4par(par))*2 max(Data(:,2)-KrLineModel4par(par))*2]);
% 
% subplot(2,2,2)
% hfit2 = plot(A.Te,(KrLineModelDiff4par(par)),...
%     'Color','Black','LineWidth',1,'LineStyle','-');
% grid on
% xlabel('Te (eV)','FontSize',14);
% ylabel('Counts/Te','FontSize',14);
% set(gca,'FontSize',14);
% set(gca,'yscale','lin');
% title(sprintf('Differential Spectrum - Doppler %s (fixed \\sigma=%.1f eV)',...
%     A.DopplerEffectFlag,A.StandDev),'FontSize',14);
% 
% subplot(2,2,4)
% [~, cor] = corrplot(errmat(1:4,1:4));
% xticklabels({'Position','Width','Amplitude','Offset'})
% yticklabels({'Position','Width','Amplitude','Offset'})
% disp(cor);
% 
% switch mypub
%     case 'ON'
%         myname = sprintf('./figures/krypton_%s_pixel%d_1.eps',...
%             TD,Pixel);
%         fprintf(2,'publish: %s\n',myname);publish_figure(fig,myname);
% end

end
