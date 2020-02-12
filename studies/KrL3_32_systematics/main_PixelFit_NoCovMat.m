addpath(genpath('../../../Samak2.0'));
clear all;
orange=1/255*[255,140,0];
Pixel = [1:148];
fitting = 'OFF';

switch fitting
    case 'OFF'
        load('FitResults_Pixels_NoCovMat.mat');
        [qU, krdatamap] = readKrData('TD',A{1,1}.TD);
        ndof = A{1,1}.nqU-4;
    case 'ON'
        A= cell(148,1);   %Kr Object for each Pixel
        parfor i= 1:148
            A{i,1} = InitKrKATRIN_krl332();
            A{i,1}.HVRipplesP2PV = 0.44;
            A{i,1}.FPD_Pixel = i;
        end
        mypar_NoCovMat = zeros(A{1,1}.nPixels,5);
        myerr_NoCovMat = zeros(A{1,1}.nPixels,5);
        mychi2min_NoCovMat = zeros(A{1,1}.nPixels,1);
        [qU, krdatamap] = readKrData('TD',A{1,1}.TD);
        
        tic  
        parfor i=1:1 %fit over all pixel  
            if max(krdatamap(i,1,:))>10   %Fit only when Bincounts are sufficient 
            sprintf('Fit of Pixel %u', i);
            [par, err, fit_xi2min] = KrL3LineFitData_NoCovMat('KrObject', A{i,1}); %Fit
            mypar_NoCovMat(i,:)= par;
            myerr_NoCovMat(i,:) = err;
            mychi2min_NoCovMat(i)= fit_xi2min; 
            else
                sprintf('No Fit of Pixel %u due to low statistic',i)
            end  
        end
        toc
     
        %Initial Values Energy and Width (same for every fit):
        LineE_i    = A{1,1}.L3_32_E_i;
        LineW_i    = A{1,1}.L3_32_W_i;
        %Initial Amplitude and Offset (different in every fit)
        LinePhi0_i_NoCovMat     = zeros(A{1,1}.nPixels,1);
        Offset_i_NoCovMat       = zeros(A{1,1}.nPixels,1);
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
        
        ndof = A{1,1}.nqU-4;
        % save results
        format long;
        clear fig; %deletes all figures
        filename_output = 'FitResults_Pixels_Simulation_NoCovMat_new.mat'
        save(filename_output); %saves whole workspace

end

%%  Display Results 

% Exclude Pixel with too low statistics for gaussian assumption

    for i = 1:148  
 if max(krdatamap(i,1,:))<10  || chi2pvalue(mychi2min_NoCovMat(i), ndof) < 0.05
    fprintf('Pixel excluded: %u \n', Pixel(i))
    E_Fit_NoCovMat(i) = NaN;
    E_FitError_NoCovMat(i) = NaN;
    W_Fit_NoCovMat(i) = NaN;
    W_FitError_NoCovMat(i) = NaN;
    Phi0_Fit_NoCovMat(i) = NaN;
    Phi0_FitError_NoCovMat(i) = NaN;
    Offset_Fit_NoCovMat(i) = NaN;
    Offset_FitError_NoCovMat(i) = NaN;;
    Pixel(i) = NaN; 
    mychi2min_NoCovMat(i) = NaN;
 end
    end   
    
% Plots 

orange=1/255*[255,140,0];
        % Plot Chi2 distribution
        f= figure(112);
        set(f, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.8]);
        set(gcf,'Color','w') %background white
        ndof = A{1,1}.nqU-4;
        plot(Pixel, mychi2min_NoCovMat./ndof,'x','Color', orange, 'MarkerSize',7);
        hold on;
        line([0 148], [1 1], 'LineStyle','--','LineWidth',2);
        hold off;
        xlabel('Pixel');
        ylabel('$\chi2/\textrm{ndof}$', 'Interpreter', 'latex');
        grid on;
        xlim([1 max(Pixel)]);
        xticks([1:9:max(Pixel)]);
        title('Pixel Fit without Covariance Matrix',...
            'FontSize', 16, 'FontWeight', 'bold');
        %export_fig './plots/WithoutCovMat/FitChi2_NoCovMat.pdf'
        %export_fig './plots/WithoutCovMat/fig/FitChi2_NoCovMat.fig'
        
        % Plot Fit Parameter Results
        f2 = figure(113);
        set(f2, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.8]);
        set(gcf,'Color','w') %background white
        subplot(2,2,1);
        errorbar(Pixel,E_Fit_NoCovMat*1e-03, E_FitError_NoCovMat*1e-03, 'x');
        xlabel('Pixel');
        ylabel('Energy (keV)');
        xlim([1 max(Pixel)]);
        xticks([0:20:max(Pixel)]);
        grid on;
        subplot(2,2,2);
        errorbar(Pixel, W_Fit_NoCovMat, W_FitError_NoCovMat, 'x');
        xlabel('Pixel');
        ylabel('Width (meV)');
        xlim([1 max(Pixel)]);
        xticks([0:20:max(Pixel)]);
        grid on;
        
        subplot(2,2,3);
        errorbar(Pixel, Phi0_Fit_NoCovMat, Phi0_FitError_NoCovMat,'x');
        xlabel('Pixel');
        ylabel('Activity (cps)');
        grid on;
        xlim([1 max(Pixel)]);
        xticks([0:20:max(Pixel)]);
        
        subplot(2,2,4);
        errorbar(Pixel, Offset_Fit_NoCovMat, Offset_FitError_NoCovMat, 'x');
        xlabel('Pixel');
        ylabel('Offset (cps)');
        grid on;
        xlim([1 max(Pixel)]);
        xticks([0:20:max(Pixel)]);
       % export_fig './plots/WithoutCovMat/FitParameter_NoCovMat.pdf'
       % export_fig './plots/WithoutCovMat/fig/FitParameter_NoCovMat.fig'
 

% Plot Spectrum

        f = figure(100);
        set(f, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.7, 0.8]);
       
        PlotPixel = 1; %Pixel of interest
        Count   = (krdatamap(PlotPixel,1,:));
        CountErr = (krdatamap(PlotPixel,2,:));
        Data = [A{1,1}.qU, Count(:), CountErr(:)];
        A{1,1}.FPD_Pixel= PlotPixel;
        A{1,1}.L3_32_Phi0_i = LinePhi0_i_NoCovMat(PlotPixel);
        A{1,1}.BKG_RateSec_i =Offset_i_NoCovMat(PlotPixel);
        A{1,1}.ComputeKrDS();
        
        figure(123);
        fig = figure('NumberTitle','off','rend','painters' ,'pos',[10 10 1300 1000]);
        set(gcf,'Color','w') %background white
        subplot(2,1,1)
        hdata = errorbar(Data(:,1)*1e-03,Data(:,2),Data(:,3),...
            's','MarkerSize',3,'Color', orange );%,'LineWidth',1);
        hold on
        hfit1 = plot(Data(:,1)*1e-03, KrLineModel4par(mypar_NoCovMat(PlotPixel,:),A{1,1}),...
            'Color','Red','LineStyle','-', 'LineWidth',1.5);
        hfit3 = line([Data(1,1)*1e-03 Data(end,1)*1e-03],[Offset_Fit_NoCovMat(PlotPixel) Offset_Fit_NoCovMat(PlotPixel)],'LineStyle','--','LineWidth',1.5);
        hold off
        grid on
        xlabel('qU (keV)','FontSize',14);
        ylabel('$\dot{\textrm{\textbf{N}}} (\textrm{\textbf{cps}})$','FontSize',14, 'Interpreter','latex');
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        %mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',...
        %   LineE_i,LineW_i);
        myfit = sprintf('Fit: \\chi2 / dof=%.1f/%.0f\n E= %.3f ± %.3f eV \n W=%.0f  ± %.0f meV \n A=%.3f  ± %.3f cps \n O=%.3f  ± %.3f cps',...
            mychi2min_NoCovMat(PlotPixel),A{1,1}.nqU-4,E_Fit_NoCovMat(PlotPixel),E_FitError_NoCovMat(PlotPixel),W_Fit_NoCovMat(PlotPixel),...
            W_FitError_NoCovMat(PlotPixel), Phi0_Fit_NoCovMat(PlotPixel),Phi0_FitError_NoCovMat(PlotPixel),Offset_Fit_NoCovMat(PlotPixel),...
            Offset_FitError_NoCovMat(PlotPixel));
        mychi2 = sprintf('Data');
        l1 = legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
        l1.FontSize = 11;
        xlim([min(A{1,1}.qU)*1e-03 max(A{1,1}.qU)*1e-03]);% 0.*min(Data(:,2)) max(Data(:,2))*1.2])
        title(sprintf('KATRIN Gaseous Krypton 83m L3-32 Pixel %u  (No Cov. Matrix)',...
           PlotPixel),'FontSize',14);
        
        subplot(2,1,2)
        plot(Data(:,1)*1e-03, (Count(:)-KrLineModel4par(mypar_NoCovMat(PlotPixel,:),A{1,1}))./CountErr(:), 'x', 'Color', orange );
        hold on
        refline = line([A{1,1}.qU(1)*1e-03, A{1,1}.qU(end)*1e-03], [0, 0], 'LineStyle', '--', 'LineWidth', 1.5)
        hold off
        grid on
        xlabel('qU (keV)','FontSize',14);
        ylabel('\textrm{\textbf{Residuals}}/$\mathbf{\sigma}$','FontSize',14, 'Interpreter','latex');
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        xlim([min(A{1,1}.qU)*1e-03 max(A{1,1}.qU)*1e-03]);%*1e-03 min(Data(:,2)-KrLineModel4par(par))*2 max(Data(:,2)-KrLineModel4par(par))*2]);
        export_fig './plots/WithoutCovMat/SpectrumPixel1_NoCovMat.pdf'
        export_fig './plots/WithoutCovMat/fig/SpectrumPixel1_NoCovMat.fig'
  
