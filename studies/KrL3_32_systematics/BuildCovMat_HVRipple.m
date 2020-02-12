%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for Building a Covariance Matrix 
% January 2018 - Lisa Schlueter MPP/TUM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MyCovMat = BuildCovMat_HVRipple(A, CovMatCompute, initFunc)

% ------ Input Parameters -----
% A = Kr Object
% CovMatCompute = Number of Spectra for CovMat
% initFunc = Function Handle(->@) for quick fit: e.g. KrL3LineFitData_NoCovMat

%Fit Results from fast Fit w/o Covariance Matrix as Initialization for Building CovMat
[quickpar, quickerr, quickfit_xi2min ndof] = initFunc('KrObject', A);
A.L3_32_E_i     = A.L3_32_E_i+quickpar(1);
A.L3_32_W_i     = A.L3_32_W_i+quickpar(2);
A.L3_32_Phi0_i  = A.L3_32_Phi0_i+quickpar(3);
A.BKG_RateSec_i = A.BKG_RateSec_i+quickpar(4);
% A.L3_32_E_i     = quickpar(1); %for hv data
% A.L3_32_W_i     = quickpar(2);
% A.L3_32_Phi0_i  = quickpar(3);
% A.BKG_RateSec_i = quickpar(4);
A.ComputeKrDS;

% Scan HVRipple Values with sigma=0.5*HVRippleP2PV
HVRipplesValue = zeros(CovMatCompute,1);
HVRipplesValue = ...     
    A.HVRipplesP2PV ...
    + 0.5*A.HVRipplesP2PV.*randn(CovMatCompute,1); 
parfor i= 1:CovMatCompute  %HV-Ripple Value cannot be <0
   if HVRipplesValue(i)<0
      HVRipplesValue(i)= 0.22 + 0.44*rand;  %uniform distributed {0.22:0.66}
   end
end
%HVRipplesValue = 2*rand(CovMatCompute,1); %uniform between 0 and 2 V
myKrIS      = zeros(CovMatCompute,A.nqU);  % for storing IS

progressbar('Generate KrIS Covariance Matrix');   
for i=1:CovMatCompute
    progressbar(i/CovMatCompute);
    A.HVRipplesP2PV = HVRipplesValue(i); 
    A.ComputeKrIS;
    myKrIS(i,:)  = model_hvdata_pixel(quickpar,A);
end
MyCovMat = cov(myKrIS);
%% Some plots for sanity check
orange=1/255*[255,140,0];   
% Plot Data + Modell with different HV-Ripples Amplitudes

    % Read Data
    if  strcmp(A.TD, 'HVL332')
        [qU krdatamap] = read_hvdata('Display','OFF','MK35',1972.4531);
    else
        [qU krdatamap] = readKrData('TD', A.TD);
    end
    Count    = (krdatamap(A.FPD_Pixel,1,:));
    CountErr = (krdatamap(A.FPD_Pixel,2,:));
    Data = [A.qU(:),Count(:),CountErr(:)];

            figS = figure('Name','Covariance Matrix Building - Spectra','NumberTitle','off','rend','painters','pos',[300 300 400 400]);
            hdataC = errorbar(Data(:,1)*1e-03,Data(:,2),Data(:,3),...
                'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold on;
            for i=1:CovMatCompute
                plot(Data(:,1)*1e-03,myKrIS(i,:),'LineWidth',0.5,'Color' ,orange);
            end
            hold off;
            grid on;
            xlabel('qU (keV)','FontSize',20);
            ylabel('$\dot{\textrm{\textbf{N}}} (\textrm{\textbf{cps}})$','FontSize',20, 'Interpreter','latex');
            %set(gca,'FontSize',14);
            set(gca,'yscale','lin');
            xlim([A.qU(1)*1e-03 A.qU(end)*1e-03]);
            dataleg = sprintf('Data KrL3-32 Pixel %u', A.FPD_Pixel);
            modelleg =sprintf('Model: HV-Ripple (0.44 ± 0.22) V');
            l1 = legend(dataleg,modelleg);
            l1.FontSize = 16;
            %title('KrL3-32 Line: Systematic Effect of HV-Ripple at KATRIN');
%export_fig './plots/WithCovMat/HVRippleVariation_Pixel1_gaussian.pdf'
%export_fig './plots/WithCovMat/fig/HVRippleVariation_Pixel1_gaussian.fig'

 % Plot Covariance Matrix
   StatErrMatrix  = (diag(CountErr(:))); %stat error only
    FitCovErrMatrix = MyCovMat + StatErrMatrix.^2;   %full CovMat
         figC = figure('Name','Covariance Matrix','NumberTitle','off','rend','painters','pos',[100 100 600 200]);
         
         subplot(1,3,1)
         imagesc(StatErrMatrix.^2);
         colormap(copper);
         colorbar;
         title('Stat. Fluct.','FontSize', 12);
         pbaspect([1 1 1])
         
         subplot(1,3,2)
         imagesc(MyCovMat);
         colormap(copper);
         colorbar;
         title('HV Ripple CM (no Stat. Fluct.)','FontSize', 12);
         pbaspect([1 1 1])
         
         subplot(1,3,3)
         imagesc(FitCovErrMatrix);
         colormap(copper);
         colorbar;
         title('HV Ripple CM + Stat. Fluct.','FontSize', 12);
         pbaspect([1 1 1])
         
%export_fig './plots/WithCovMat/CovarianceMatrix_Pixel1_gaussian.pdf'
%export_fig './plots/WithCovMat/fig/CovarianceMatrix_Pixel1_gaussian.fig'

end