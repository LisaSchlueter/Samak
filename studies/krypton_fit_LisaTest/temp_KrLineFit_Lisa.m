function [par err chi2min ndof] = KrLineFit_Lisa(varargin)
%
%            FIT INTEGRAL SPECTRUM
%              Fit a 83mKr Lines
%
%          Th. Lasserre - CEA Saclay
%                August 2017
%

% Initialization
clear par ; clear parnomix;
addpath(genpath('../../../Samak2.0'));
% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub','YES',@(x)ismember(x,{'ON','OFF','YES'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Mode','MoSFitterJuly2017',@(x)ismember(x,{'Data','Sim','MoSFitterJuly2017'}));
p.addParameter('CPS','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('FPD_Segmentation','OFF',@(x)ismember(x,{'OFF','RING','PIXEL'}));
p.addParameter('Pixel',1,@(x)isfloat(x) && x>0);
p.addParameter('nfit',1,@(x)isfloat(x) && x>=0);
p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
p.addParameter('CovMat','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TD','MoSFitterJuly2017',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn','Flat20','Flat30','Flat50','Flat60','Flat100','KrK32','KrL3_32','KrL3_32_HS','MoSFitterJuly2017','KrL3_32_Satellites'}));
p.addParameter('DopplerEffectFlag','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('SystCorr','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('HVRipples','OFF',@(x)ismember(x,{'OFF','ON'}));
p.addParameter('HVRipplesP2PV',0.52,@(x)isfloat(x) && x>0); % V or eV
p.addParameter('ConvFlag','OFF', @(x)ismember(x,{'DopplerEffect','Voigt','OFF'}))
p.addParameter('MultiPeaksFlag','OFF',@(x)ismember(x,{'OFF','ON'}));

p.parse(varargin{:});
fign       =    p.Results.fign;
pub        =    p.Results.pub;
display    =    p.Results.display;
Mode       =    p.Results.Mode;
CPS        =    p.Results.CPS;
FPD_Segmentation  = p.Results.FPD_Segmentation;
Pixel      =    p.Results.Pixel;
nfit       =    p.Results.nfit;
Fitter     =    p.Results.Fitter;
CovMat     =    p.Results.CovMat;
TD         =    p.Results.TD;
DopplerEffectFlag = p.Results.DopplerEffectFlag;
SystCorr   =    p.Results.SystCorr;
HVRipples  = p.Results.HVRipples;
HVRipplesP2PV = p.Results.HVRipplesP2PV;
ConvFlag   = p.Results.ConvFlag;
MultiPeaksFlag = p.Results.MultiPeaksFlag;

clear FitCovErrMatrix;
global FitCovErrMatrix;

% Parametrization: True Value
global A ;

switch TD
    case {'KrK32','KrL3_32','KrL3_32_HS'}
        A=InitKrKATRIN_LisaFit('TD',TD,'FPD_Pixel',Pixel,...
            'ConvFlag',ConvFlag,...
            'CPS',CPS,'FPD_Segmentation',FPD_Segmentation,...
            'HVRipples',HVRipples,'HVRipplesP2PV',HVRipplesP2PV);
    case 'MoSFitterJuly2017'
        A=InitKrKATRIN_LisaFit('TD',TD,'FPD_Segmentation','OFF',...
            'HVRipples','ON',...
            'ConvFlag',ConvFlag,'CPS','OFF');
     case  'KrL3_32_Satellites'
            A=InitKrKATRIN_LisaFit('TD',TD,'FPD_Segmentation','OFF',...
                'HVRipples','ON','ConvFlag',ConvFlag,'CPS','OFF', ...
                'MultiPeaksFlag', MultiPeaksFlag, 'DopplerEffectFlag', DopplerEffectFlag);
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
   case 'KrL3_32_Satellites'
        LineE_i    = A.L3_32_E_i;
        LineW_i    = A.L3_32_W_i;
        LinePhi0_i = A.L3_32_Phi0_i;
        LineE_s3_i = A.L3_32_E_s3_i; %Theoretical/Initial Energy Satellite Line
        LinePhi0_s3_i = A.L3_32_Phi0_s3_i;  %Intensity Satellite line (relative)
end

switch CovMat
    case 'ON'
        nCMTrials = 10000;
        fName = sprintf('../../data/CovMat/WGTSMACE_CovMat_%gTrials_%s_Year_%g.mat',...
            nCMTrials,A.TD,A.TimeYear);
        if exist(fName, 'file') == 2
            CM=importdata(fName);
        else
            %A.ComputeTBDISCovarianceMatrix('nTrials',nCMTrials);
            CM=importdata(fName);
        end
end

% Loop on fits
switch MultiPeaksFlag
        case 'ON'
           npar=7;  
           i_E_s3 = 0;  %Initialization satellite
           i_Phi0_s3 = 0; 
        case 'OFF'
           npar       = 5;
    end 

fit_p      = ones(npar,nfit);
fit_xi2min = ones(1,nfit);

% Init
i_E        = 0; i_W        = 0;
i_Bkg      = 0; i_Phi0     = 0;
i_N        = 0;

switch Fitter
    case 'Minuit'
    case 'Matlab'
        %             tbdisf = @(mb,qb,bb,nb,qu) A.ComputeTBDISf(qu,mb,qb,bb,nb,0,0);
        %             tmp = @(mb,qb,bb,nb,e) interp1(A.qU,tbdisf(mb,qb,bb,nb,A.qU),e);
        %             myf = fittype(@(mb,qb,bb,nb,qu) tmp(mb,qb,bb,nb,qu),...
        %                 'independent','qu','coefficients', {'mb','qb','bb','nb'});
        %             opts = fitoptions('StartPoint',[0,0,0,0],...
        %                 'Method', 'NonlinearLeastSquares','Display','off',...
        %                 'Weights',1./A.TBDIS);
end

% Data / Sim
switch Mode
    case 'Data'
        [qU krdatamap] = readKrData('TD',TD);
    case 'MoSFitterJuly2017'
        [qU krdatamap] = readKrData('TD',TD);
end

% Background / Amplitude Initialization
switch TD
    case 'KrK32'
    case {'KrL3_32'}
        qumin=1;
        tmp = krdatamap(Pixel,1,qumin);
        A.L3_32_Phi0_i  =tmp;
        A.BKG_RateSec_i = krdatamap(Pixel,1,end);
    case {'KrL3_32_HS'}
        qumin=100;
        tmp = krdatamap(Pixel,1,qumin);
        A.L3_32_Phi0_i  = tmp;
        A.BKG_RateSec_i = krdatamap(Pixel,1,end);
    case 'MoSFitterJuly2017'
        tmp = 0;
     case 'KrL3_32_Satellites'
            qumin=1;
            tmp = krdatamap(Pixel,1,qumin);
            A.L3_32_Phi0_i  = tmp;
            A.BKG_RateSec_i = krdatamap(Pixel,1,end);    
end

% Number of Fits
tic
progressbar('Krypton Line Fits');
for f=1:1:nfit
    progressbar(f/nfit);
    
    switch Mode
        case 'Data'
            Count    = (krdatamap(Pixel,1,qumin:end));
            CountErr = (krdatamap(Pixel,2,qumin:end));
            Data = [A.qU(:),Count(:),CountErr(:)];
            %disp(Data);
        case 'Sim'
            A.ComputeKrDS(); A.ComputeKrIS();
            AddStatFluctKrIS(A);
            %AddStatSystFluctTBDIS(A,'CovMat',CovMat,'SystCorr',SystCorr);
            Data = [A.qU,A.KrIS,A.KrISE];
        case 'MoSFitterJuly2017'
            Count    = (krdatamap(1,1,:));
            CountErr = (krdatamap(1,2,:));
            Data = [A.qU(:),Count(:),CountErr(:)];
            %disp(Data);
    end
    
    switch Fitter
        case 'Minuit'
            % Initializing Fit
          switch MultiPeaksFlag
              case 'ON'
                  ParIni = [i_E i_W i_Phi0 i_Bkg i_N i_E_s3 i_Phi0_s3];
                  parnames = ['E W Phi0 Bkg N E_s3 Phi0_s3'];
                        
              case 'OFF'
                  ParIni = [i_E i_W i_Phi0 i_Bkg i_N];
                   parnames = ['E W Phi0 Bkg N'];
          end 
          
            tmparg = sprintf(['set pri -10 ; fix 2 3 4 5'...
                'set now; min ; imp' ]);
            Args = {ParIni, Data, '-c',tmparg};
            
            switch CovMat
                case 'OFF'
                    fprintf(2,'No Systematics - Stat. + Pulls\n');
                  switch 'MultiPeaksFlag'
                           case 'OFF'
                        [par, err, chi2min, errmat] = fminuit('KrLineChi2',Args{:});
                            case 'ON' 
                         [par, err, chi2min, errmat] = fminuit('KrLineChi2_sat',Args{:}); 
                  end 
                % ndof = size(A.qU)-npar; %new Lisa
                case 'ON'
                    fprintf(2,'Stat + Systematics (Covariance Matrix)\n');
                    switch SystCorr
                        case 'OFF'
                            FitCovErrMatrix = diag(diag(CM.WGTSMACE_CovMat)) + diag(A.KrIS);
                        case 'ON'
                            FitCovErrMatrix = (CM.WGTSMACE_CovMat) + diag(A.KrIS);
                    end
                    %[ par, err, chi2min, errmat ] = fminuit('KrLineChi2CovMat',Args{:});
            end
            
           fit_p(:,f)=par;  fit_xi2min(f) = chi2min;
            
        case 'Matlab'
            %                 [fit1,gof,fitinfo] = fit(A.qU,A.TBDIS,myf,opts);
            %                 par = [fit1.mb fit1.qb fit1.bb fit1.nb 0 0];
            %                 fit_p(:,f)=par;
            %                 err = [0 0 0 0 0 0];
            %                 chi2min = gof.sse; fit_xi2min(f) = chi2min;
    end
    
    switch display
        case 'ON'
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf('  Processing \t= %g %% \n',f/nfit*100);
            fprintf(2,'--------------------------------------------------------------\n');
            fprintf(2,'  E \t= %.3f \t� \t %g \t eV \n', (LineE_i+par(1)),err(1));
            fprintf(2,'  W \t= %.3f \t� \t %g \t meV \n', (LineW_i+par(2))*1e3,err(2)*1e3);
            %fprintf(2,'  Phi0 \t= %.3f \t� \t %g \t    \n', (tmp+par(3)),err(3));
            fprintf(2,'  Bkg \t= %.3f \t� \t %g \t cps \n', (A.BKG_RateSec_i+par(4)),err(4));
            fprintf(2,'  N \t= %g \t� \t %g \t    \n', par(5),err(5));
            ndof = A.nqU-4;
            fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
            fprintf(2,'--------------------------------------------------------------\n');
            
    end
end
toc
end
% %% Plot Results
% figure(fign+999)
% subplot(2,1,1)
% hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
%     'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
% hold on
% %     hfit1 = plot(A.Te,KrLineModelDiff4par(par)./trapz(A.Te,KrLineModelDiff4par(par)).*trapz(A.qU,KrLineModel4par(par))/4,...
% %         'LineWidth',1,'LineStyle','--','Color','Black');
% hfit1 = plot(Data(:,1),KrLineModel4par(par),...
%     'Color','Black','LineWidth',1,'LineStyle','-');
% hold off
% grid on
% xlabel('qU (eV)','FontSize',10);
% ylabel('Counts','FontSize',10);
% switch Mode
%     case 'MoSFitterJuly2017'
%         title(sprintf('KATRIN Krypton - %s - Doppler %s',A.TD,A.DopplerEffectFlag));
%     case {'DATA','Sim'}
%         title(sprintf('KATRIN Krypton - %s - %s - Pixel %.0f - Doppler %s',A.TD,Mode,Pixel,A.DopplerEffectFlag));
% end
% set(gca,'FontSize',12);
% set(gca,'yscale','lin');
% mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',LineE_i,LineW_i);
% myfit = sprintf('Fit: E=%.3f \\pm %.3f eV - W=%.3f \\pm %.3f eV',LineE_i+par(1),err(1),LineW_i+par(2),err(2));
% mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-4);
% legend([hdata  hfit1],mychi2,myfit,'Location','NorthEast') ;% legend(a,'boxoff');
% axis([min(A.qU) max(A.qU)+1 0.7*min(Data(:,2)) max(Data(:,2))*1.2])
% 
% subplot(2,1,2)
% hdata = errorbar(Data(:,1),Data(:,2)-KrLineModel4par(par),Data(:,3),...
%     'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
% hold on
% line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--');
% hold off
% grid on
% xlabel('qU (eV)','FontSize',10);
% ylabel('Residuals','FontSize',10);
% %title(sprintf('KATRIN Integral - %s - %s - Pixel %.0f',A.TD,Mode,Pixel));
% set(gca,'FontSize',12);
% set(gca,'yscale','lin');
% %     mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',LineE_i,LineW_i);
% %     myfit = sprintf('Fit: E=%.2f \\pm %.2f eV - W=%.2f \\pm %.2f eV',LineE_i+par(1),err(1),LineW_i+par(2),err(2));
% %     mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-4);
% %     legend([hdata  hfit1],mychi2,myfit,'Location','NorthEast') ;% legend(a,'boxoff');
% %axis([min(A.qU) max(A.qU)+1 0.7*min(Data(:,2)) max(Data(:,2))*1.2])
% axis([min(A.qU) max(A.qU)+1 min(Data(:,2)-KrLineModel4par(par))*2 max(Data(:,2)-KrLineModel4par(par))*2])
% 
% pub = 'YES';
% switch pub
%     case 'YES'
%         switch Mode
%             case 'MoSFitterJuly2017'
%                 myname = sprintf('./figures/MoSFitterJuly2017_1.eps');
%                 fprintf(2,'publish: %s\n',myname);publish_figure(999+fign,myname);
%             case {'DATA','Sim'}
%                 myname = sprintf('./figures/kryptonHS_%s_pixel%d_1.eps',...
%                     TD,Pixel);
%                 fprintf(2,'publish: %s\n',myname);publish_figure(999+fign,myname);
%         end
%         
%         
%         %             figure(fign+1000)
%         %             [h cor] = corplot(errmat(1:4,1:4));
%         %             xticklabels({'Position','Width','Amplitude','Background'})
%         %             yticklabels({'Position','Width','Amplitude','Background'})
%         %             title(sprintf('KATRIN Krypton - Correlation Matrix - %s - Doppler %s\n',A.TD,A.DopplerEffectFlag));
%         %             fprintf(2,'Error Matrix�\n');
%         %             disp(cor);
%         
%         switch Mode
%             case 'MoSFitterJuly2017'
%                 switch pub
%                     case 'YES'
%                         myname = sprintf('./figures/MoSFitterJuly2017_2.eps');
%                         fprintf(2,'publish: %s\n',myname); %publish_figure(1000+fign,myname);
%                 end
%             case {'DATA','Sim'}
%                 myname = sprintf('./figures/krypton_%s_pixel%d_1.eps',...
%                     TD,Pixel);
%                 fprintf(2,'publish: %s\n',myname);publish_figure(999+fign,myname);
%         end
%         
%         if nfit<5
%             % For Saving:
%             par(1) = LineE_i         + par(1);
%             par(2) = LineW_i         + par(2);
%             par(3) = tmp             + par(3);
%             par(4) = A.BKG_RateSec_i + par(4);
%         end
%         
%        % return;
% 
% 
%     figure(fign)
%     subplot(2,1,1)
%     hfit1 = plot(Data(:,1),KrLineModel4par(par),...
%         'LineWidth',1,'LineStyle','-','Color','Black');
%     hold on;
%     hfit2 = plot(Data(:,1),KrLineModel4par(par),...
%         'LineWidth',1,'LineStyle','-','Color','Black');
%     hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
%         'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
%     %errorbar_tick(hdata,200);
%     hold off
%     grid on
%     xlabel('qU (eV)','FontSize',10);
%     ylabel('Counts','FontSize',10);
%     title(sprintf('KATRIN Integral - %s - %g y - %g mcps',A.TD,A.TimeYear,A.BKG_RateSec_i));
%     set(gca,'FontSize',12);
%     set(gca,'yscale','lin');
%     mydata = sprintf('Data: E=%.2f eV \n',LineE_i);
%     myfit = sprintf('Fit: E=%.2f \\pm %.2f eV',LineE_i+par(1),err(1));
%     mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU-4);
%     legend([hdata hfit1 hfit2],mydata,myfit,mychi2,'Location','NorthEast') ;% legend(a,'boxoff');
%     axis([min(A.qU) max(A.qU)+1 min(Data(:,2)) max(Data(:,2))*1.1])
%     
%     subplot(2,1,2)
%     parnomix = zeros(1,4,1); parnomix = par; parnomix(1) = 0; parnomix(2) = 0;
%     hfit = plot(Data(:,1),...
%         KrLineModel4par(par)./KrLineModel4par(parnomix)-1,...
%         'Color','Black','LineWidth',1,'LineStyle','-');
%     hold on;
%     hdata = errorbar(Data(:,1),...
%         Data(:,2)./KrLineModel4par(parnomix)-1,...
%         Data(:,3)./KrLineModel4par(parnomix),...
%         'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
%     %errorbar_tick(hdata,200);
%     hold off;
%     grid on
%     xlabel('qU (eV)','FontSize',10);
%     ylabel('Spectral distorsion','FontSize',10);
%     set(gca,'FontSize',12);
%     a= legend([hdata hfit],'Data/No Mixing - 1','Fit/No Mixing - 1','Location','NorthWest');
%     %legend(a,'boxoff');
%     %axis([min(A.qU) max(A.qU)+1 min((Data(:,2)./KrLineModel4par(parnomix)-1)) max(Data(:,3)./KrLineModel4par(parnomix))*3])
% 
%     switch pub
%         case 'YES'
%            myname = sprintf('./figures/f1-krkatrin_%s_pixel%.0f.eps',...
%                 TD,Pixel); 
%             fprintf(2,'publish: %s\n',myname);publish_figure(fign,myname);
%     end
%     
%     figure(fign+1)
%     subplot(2,2,1);
%     title('K32_E','FontSize',10)
%     nhist(fit_p(1,:)*1e3,'text','pdf','color','sequential');xlabel('meV','FontSize',8); grid on
%     subplot(2,2,2);
%     title('K32_W','FontSize',12)
%     nhist((fit_p(2,:))*1e3,'text','pdf','color','sequential');xlabel('meV','FontSize',8); grid on
%     subplot(2,2,3);
%     title('Normalzation','FontSize',12)
%     nhist(fit_p(3,:),'text','pdf','color','sequential');xlabel('no unit','FontSize',8); grid on
%     subplot(2,2,4);
%     title('Background','FontSize',12)
%     nhist(((fit_p(4,:))),'text','pdf','color','sequential');xlabel('mcps','FontSize',8); grid on
%     switch pub
%         case 'ON'
%            myname = sprintf('./figures/f2-krkatrin_%g-mcps_%g-numsq.eps',...
%                 A.BKG_RateAllFPDSec*1e3,0); 
%             fprintf(2,'publish: %s\n',myname);publish_figure(fign+1,myname);
%     end
%     
%     figure(fign+2)
%     ndhist(fit_p(1,:)*1e3,((fit_p(2,:)))*1e3);
%     colorbar
%     ylabel('E (meV)','FontSize',10);
%     xlabel('W (meV)','FontSize',10);
%     R1 = corrcoef(fit_p(1,:),((fit_p(2,:))));
%     switch pub
%         case 'ON'
%            myname = sprintf('./figures/f3-krkatrin_%g-mcps_%g-numsq.eps',...
%                 A.BKG_RateAllFPDSec*1e3,0); 
%             fprintf(2,'publish: %s\n',myname);publish_figure(fign+2,myname);
%     end
%     
%     
%     % Display
%     % A.DisplayKrInfo();
%     
%                 
% end
