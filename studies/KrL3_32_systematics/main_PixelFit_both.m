%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fit of all Pixel with and without Covariance Matrix
% 
% Lisa Schlueter (January 2018) MPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
format long; 
addpath(genpath('../../../Samak2.0'));
orange=1/255*[255,140,0];        
  
fittingCovMat = 'OFF';
fittingNoCovMat = 'ON';

switch fittingCovMat
    %case 'OFF'
       % load('FitResults_Pixels_CovMat.mat');
        %[qU krdatamap] = readKrData('TD',A{1,1}.TD);         
    case 'ON' 
        %Kr Object for each Pixel
        A= cell(148,1); %CovMat Object
        parfor i= 1:148
            A{i,1} = InitKrKATRIN_krl332();
            A{i,1}.HVRipplesP2PV = 0.44;   %central value
            A{i,1}.FPD_Pixel = i;
        end
        mypar = zeros(A{1,1}.nPixels,5); %Init Fit Results
        myerr = zeros(A{1,1}.nPixels,5);
        mychi2min = zeros(A{1,1}.nPixels,1);
        [qU krdatamap] = readKrData('TD',A{1,1}.TD);
        
        tic 
        parfor i = 1:148   
            if max(krdatamap(i,1,:))>10   %Fit only when Bincounts are sufficient
                sprintf('Fit of Pixel %u', i);
                [par, err, chi2min, ndof] = KrL3LineFitData_CovMat('CovMat','ON',...
                    'CovMatCompute',1e4, 'KrObject', A{i,1});
                mypar(i,:)= par;
                myerr(i,:) = err;
                mychi2min(i)= chi2min;
            else
                sprintf('No Fit of Pixel %u due to low statistic',i)
            end
        end
        toc
  %Save Fit Results   
     %Initial Values Energy and Width (NOT!!!! same for every fit):
        LineE_i    = A{1,1}.L3_32_E_i;
        LineW_i    = A{1,1}.L3_32_W_i;
     %Initial Amplitude and Offset (different in every fit)
        LinePhi0_i     = zeros(A{1,1}.nPixels,1);
        Offset_i       = zeros(A{1,1}.nPixels,1);
        parfor i=1:148
            LinePhi0_i(i) = krdatamap(i,1,1)-krdatamap(i,1,end);
            Offset_i(i)   = krdatamap(i,1,end);
        end
 
        E_Fit             = (LineE_i+mypar(:,1));
        E_FitError        = (myerr(:,1));
        
        W_Fit             = (LineW_i+mypar(:,2))*1e3;
        W_FitError        = (myerr(:,2))*1e3;;
        
        Phi0_Fit          = (LinePhi0_i(:)+mypar(:,3));
        Phi0_FitError     = myerr(:,3);
        
        Offset_Fit        = (Offset_i+mypar(:,4));
        Offset_FitError   = myerr(:,4);
        
        %Save results
        Pixel = [1:148]; 
        clear fig; %deletes all figures
        filename_outputCovMat = 'FitResults_Pixels_CovMat_new.mat'
        save(filename_outputCovMat);%, E_Fit, E_FitError, W_Fit,W_FitError,Phi0_Fit, Phi0_FitError,...
            %Offset_Fit,Offset_FitError, LinePhi0_i, Offset_i,LineE_i, LineW_i,...
            %mychi2min, mypar, myerr, A); %save fit results   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit Without CovMat %%%%%%%%%%%%%%%%%%%%%%%%%%
switch fittingNoCovMat
    case 'OFF'
        load('FitResults_Pixels_NoCovMat.mat');
        [qU, krdatamap] = readKrData('TD',B{1,1}.TD);
        ndof = B{1,1}.nqU-4;
    case 'ON'
        B= cell(148,1);   %Kr Object for each Pixel
        parfor i= 1:148
            B{i,1} = InitKrKATRIN_krl332();
            B{i,1}.HVRipplesP2PV = 0.44;
            B{i,1}.FPD_Pixel = i;
        end
        mypar_NoCovMat = zeros(B{1,1}.nPixels,5);
        myerr_NoCovMat = zeros(B{1,1}.nPixels,5);
        mychi2min_NoCovMat = zeros(B{1,1}.nPixels,1);
        [qU, krdatamap] = readKrData('TD',B{1,1}.TD);
        
        tic  
        parfor i=1:148 %fit over all pixel  
            if max(krdatamap(i,1,:))>10   %Fit only when Bincounts are sufficient 
            sprintf('Fit of Pixel %u (No CovMat)', i);
            [par, err, fit_xi2min] = KrL3LineFitData_NoCovMat('KrObject', B{i,1}); %Fit
            mypar_NoCovMat(i,:)= par;
            myerr_NoCovMat(i,:) = err;
            mychi2min_NoCovMat(i)= fit_xi2min; 
            else
                sprintf('No Fit of Pixel %u due to low statistic',i)
            end  
        end
        toc
     
        %Initial Values Energy and Width (same for every fit):
        LineE_i    = B{1,1}.L3_32_E_i;
        LineW_i    = B{1,1}.L3_32_W_i;
        %Initial Amplitude and Offset (different in every fit)
        LinePhi0_i_NoCovMat     = zeros(B{1,1}.nPixels,1);
        Offset_i_NoCovMat       = zeros(B{1,1}.nPixels,1);
        parfor i=1:148
            LinePhi0_i_NoCovMat(i) = krdatamap(i,1,1)-krdatamap(i,1,end);
            Offset_i_NoCovMat(i)   = krdatamap(i,1,end);
        end
        
        %Fit Results: Intial Value + Fit Parameter
        E_Fit_NoCovMat       = (LineE_i+mypar_NoCovMat(:,1));
        E_FitError_NoCovMat  = (myerr_NoCovMat(:,1));
        
        W_Fit_NoCovMat       = (LineW_i+mypar_NoCovMat(:,2))*1e3;
        W_FitError_NoCovMat  = (myerr_NoCovMat(:,2))*1e3;
        
        Phi0_Fit_NoCovMat       = (LinePhi0_i_NoCovMat(:)+mypar_NoCovMat(:,3));
        Phi0_FitError_NoCovMat  = myerr_NoCovMat(:,3);
        
        Offset_Fit_NoCovMat      = (Offset_i_NoCovMat(:)+mypar_NoCovMat(:,4));
        Offset_FitError_NoCovMat = myerr_NoCovMat(:,4);
        
        ndof = B{1,1}.nqU-4;
        % save results
        format long;
        clear fig; %deletes all figures
        filename_output = 'FitResults_Pixels_NoCovMat_new.mat'
        save(filename_output)%, E_Fit_NoCovMat, E_FitError_NoCovMat, W_Fit_NoCovMat,W_FitError_NoCovMat,...
            %Phi0_Fit_NoCovMat, Phi0_FitError_NoCovMat,Offset_Fit_NoCovMat,Offset_FitError_NoCovMat,...
            %LinePhi0_i_NoCovMat, Offset_i_NoCovMat,LineE_i, LineW_i,...
            %mychi2min_NoCovMat, ndof, mypar_NoCovMat, myerr_NoCovMat, B); %save fit results 
end

% Exclude Pixel with too low statistics
% Pixel = [1:148];  
%     parfor i= 1:148 
%  if max(krdatamap(i,1,:))<10  || chi2pvalue(mychi2min(i), ndof) < 0.05 || chi2pvalue(mychi2min_NoCovMat(i), ndof) < 0.05
%     fprintf('Pixel excluded: %u \n', Pixel(i))
%     E_Fit(i) = NaN;
%     E_FitError(i) = NaN;
%     W_Fit(i) = NaN;
%     W_FitError(i) = NaN;
%     Phi0_Fit(i) = NaN;
%     Phi0_FitError(i) = NaN;
%     Offset_Fit(i) = NaN;
%     Offset_FitError(i) = NaN;
%     Pixel(i) = NaN; 
%     mychi2min(i) = NaN; 
%     
%     E_Fit_NoCovMat(i) = NaN;
%     E_FitError_NoCovMat(i) = NaN;
%     W_Fit_NoCovMat(i) = NaN;
%     W_FitError_NoCovMat(i) = NaN;
%     Phi0_Fit_NoCovMat(i) = NaN;
%     Phi0_FitError_NoCovMat(i) = NaN;
%     Offset_Fit_NoCovMat(i) = NaN;
%     Offset_FitError_NoCovMat(i) = NaN;
%     mychi2min_NoCovMat(i) = NaN;
%  end
%     end 
%     
%% PLOTS 
% f = figure(100); %Plot Chi2 distribution
% set(f, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.8]);  
% set(gcf,'Color','w') %background white
% plot(Pixel, mychi2min./ndof,'x', 'Color', orange, 'MarkerSize',7); 
% hold on;
% line([0 148], [1 1], 'LineStyle','--','LineWidth',2);
% hold off;
% xlabel('Pixel', 'FontSize', 16);
% ylabel('$\mathbf{\chi2}/\textrm{\textbf{ndof}}$', 'Interpreter', 'latex', 'FontSize', 16);
% grid on;
% xlim([1 max(Pixel)]);
% xticks([1:9:max(Pixel)]);
% title('Pixel Fit with Covariance Matrix',...
%     'FontSize', 16, 'FontWeight', 'bold');
% %export_fig './plots/WithCovMat/FitChi2_CovMat.pdf'
% %export_fig './plots/WithCovMat/fig/FitChi2_CovMat.fig'
% 
% %FPDViewer(mychi2min./ndof);
% %FPDViewer(E_Fit);   
% %FPDViewer(mychi2min./ndof);
% %export_fig './plots/WithCovMat/FPD_WError.png'
% %export_fig './plots/WithCovMat/FPD_chi2min_ndof.png'
% %% Plot Spectrum
% PlotPixel = 1; %Pixel of interest
% Count   = (krdatamap(PlotPixel,1,:));
% CountErr = (krdatamap(PlotPixel,2,:));
% Data = [qU , Count(:), CountErr(:)];
% A{PlotPixel,1}.FPD_Pixel= PlotPixel;
% A{PlotPixel,1}.L3_32_Phi0_i = LinePhi0_i(PlotPixel);
% A{PlotPixel,1}.BKG_RateSec_i =Offset_i(PlotPixel);
% A{PlotPixel,1}.ComputeKrDS;
% 
% fig = figure('NumberTitle','off','rend','painters' ,'pos',[10 10 1300 1000]);
% set(gcf,'Color','w') %background white
% subplot(2,1,1)
% hdata = errorbar(Data(:,1)*1e-03,Data(:,2),Data(:,3),...
%     's','MarkerSize',3,'Color', orange );%,'LineWidth',1);
% hold on
% hfit1 = plot(Data(:,1)*1e-03, krl332minuit_modelint(mypar(PlotPixel,:), A{PlotPixel,1} ),...
%     'Color','Red','LineStyle','-', 'LineWidth',1.5);
% hfit3 = line([Data(1,1)*1e-03 Data(end,1)*1e-03],[Offset_Fit(PlotPixel) Offset_Fit(PlotPixel)],'LineStyle','--','LineWidth',1.5);
% hold off
% grid on
% xlabel('qU (keV)','FontSize',14);
% ylabel('$\dot{\textrm{\textbf{N}}} (\textrm{\textbf{cps}})$','FontSize',14, 'Interpreter','latex');
% set(gca,'FontSize',14);
% set(gca,'yscale','lin');
% %mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
% %   LineE_i,LineW_i);
% myfit = sprintf('Fit: \\chi2 / dof=%.1f/%.0f\n E= %.3f ± %.3f eV \n W=%.0f  ± %.0f meV \n A=%.3f  ± %.3f cps \n O=%.3f  ± %.3f cps',...
%     mychi2min(PlotPixel),ndof,E_Fit(PlotPixel),E_FitError(PlotPixel),W_Fit(PlotPixel),...
%     W_FitError(PlotPixel), Phi0_Fit(PlotPixel),Phi0_FitError(PlotPixel),Offset_Fit(PlotPixel),...
%     Offset_FitError(PlotPixel));
% mychi2 = sprintf('Data');
% l1 = legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
% l1.FontSize = 11;
% xlim([min(A{1,1}.qU)*1e-03 max(A{1,1}.qU)*1e-03]);% 0.*min(Data(:,2)) max(Data(:,2))*1.2])
% title(sprintf('KATRIN Gaseous Krypton 83m -{%s} Pixel %u Data and Fit (Cov. Matrix)',...
%     'KrL3-32', PlotPixel),'FontSize',14);
% 
% subplot(2,1,2)
% plot(Data(:,1)*1e-03, (Count(:)-krl332minuit_modelint(mypar(PlotPixel,:),A{1,1}))./CountErr(:), 'x', 'Color', orange );
% hold on
% refline = line([A{1,1}.qU(1)*1e-03, A{1,1}.qU(end)*1e-03], [0, 0], 'LineStyle', '--', 'LineWidth', 1.5);
% hold off
% grid on
% xlabel('qU (keV)','FontSize',14);
% ylabel('\textrm{\textbf{Residuals}}/$\mathbf{\sigma}$','FontSize',14, 'Interpreter','latex');
% set(gca,'FontSize',14);
% set(gca,'yscale','lin');
% xlim([min(A{1,1}.qU)*1e-03 max(A{1,1}.qU)*1e-03]);%*1e-03 min(Data(:,2)-KrLineModel4par(par))*2 max(Data(:,2)-KrLineModel4par(par))*2]);
% export_fig './plots/WithCovMat/fig/SpectrumPixel33_CovMat.fig'
% export_fig './plots/WithCovMat/SpectrumPixel33_CovMat.pdf'