addpath(genpath('../../../Samak2.0'));

% Settings
saveplot = 'OFF';
chi2 = 'Stat'; %CM, CMFrac, Stat
VFTRuns = [40257,40263];%, 40263,40264];%:40260,40263:40266];
fix = 0;
% Excluding Data Points?
exclData = 'ON';
exclnPoints = 12; %exclude up to (exclnPoints-1) Points
exclDataStop = 9999; % 9999== take all point until qUmin
if strcmp(exclData,'OFF')
    exclnPoints = 1;
end
% Column Density Scan?
rhodScan = 'OFF';

for r =1:numel(VFTRuns)
RunNr = VFTRuns(r);
for i=1:exclnPoints 
exclDataStart = i;

%% Data 
D = importdata([num2str(RunNr),'.mat']);
switch RunNr
    case {40259,40260,40263, 40264, 40265,40266}
        Data = [D.qU, D.TBDIS, sqrt(D.TBDIS)];
    case {40257,40258}
        Data = [D.qU(3:end), D.TBDIS(3:end), sqrt(D.TBDIS(3:end))];
end
%% Model
switch rhodScan
    case 'ON'
        WGTS_CD_MolPerCm2_local = D.WGTS_CD_MolPerCm2.*[0.2:0.1:1.5];
        parRD = zeros(4,numel(WGTS_CD_MolPerCm2_local));
        errRD = zeros(4,numel(WGTS_CD_MolPerCm2_local));
        chi2RD = zeros(numel(WGTS_CD_MolPerCm2_local),1);
    case 'OFF'
        WGTS_CD_MolPerCm2_local = D.WGTS_CD_MolPerCm2;
end
for rd =1:1:numel(WGTS_CD_MolPerCm2_local) %rhod scan
M = ref_runsummaries(RunNr,'ISCS','Theory',...
                    'recomputeRF','OFF','WGTS_CD_MolPerCm2',WGTS_CD_MolPerCm2_local(rd)); 
M.ComputeTBDDS; M.ComputeTBDIS;

%% Covariance Matrix
cm = importdata('CM_VFT2ud_20percent2.mat');
%cm = importdata('CM_VFT2ud.mat');
CM_Stat = diag(Data(:,2));
CMFrac_Stat = diag(1./Data(:,2));
FitCM = Data(:,2).*cm.covmatfrac.*Data(:,2)' + CM_Stat;
FitCMFrac = cm.covmatfrac + CMFrac_Stat; 

%% Fit
tic;
switch chi2
    case 'CM'
        F = FITC('SO',M,'DATA',Data,'fitter','minuit',...
            'chi2name','chi2CM', 'COVMAT', FitCM,...
            'pulls',[4 inf inf inf],...
            'exclDataStart',exclDataStart); % excludes start-1 points
    case 'CMFrac'
        F = FITC('SO',M,'DATA',Data,'fitter','minuit',...
            'chi2name','chi2CMFrac', 'COVMATFrac', FitCMFrac,...
            'pulls',[4 inf inf inf],...
            'exclDataStart',exclDataStart);
    case 'Stat'
        F = FITC('SO',M,'DATA',Data,'fitter','minuit',...
            'chi2name','chi2',...
            'pulls',[4 inf inf inf],...
            'exclDataStart',exclDataStart);
        FitCM = CM_Stat;
end
 toc;

par = F.RESULTS{1};
err = F.RESULTS{2};
chi2min = F.RESULTS{3};
dof = numel(M.qU(exclDataStart:end))-4+fix;
switch rhodScan
    case 'ON'
      parRD(:,rd) =  par;
      errRD(:,rd) =  err;
      chi2RD(rd)  = chi2min;      
end
%% Plots Spectrum + Fit
fig1 = figure(5+rd);
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,1]);
subplot(2,1,1);
hfit1 = boundedline(M.qU,M.TBDIS,diag(sqrt(FitCM)),'alpha');
hold on;
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
%hfit2  = plot(M.qU,M.TBDIS,...
 %   'Color','Black','LineWidth',1,'LineStyle','-');
%hline = line([Data(1,1) Data(end,1)],[(M.BKG_RateSec_i+par(3))*M.TimeSec*M.qUfrac(end) (M.BKG_RateSec_i+par(3))*M.TimeSec*M.qUfrac(end)],'LineStyle','--','Color','Red');
hold off;
set(gca, 'YScale', 'log')
xlabel('qU (eV)');
ylabel('counts')
xlim([min(M.qU(exclDataStart:end)) max(M.qU)]);
myfitleg = sprintf(' Fit: \\chi2 / dof=%.2f/%.0f \n m^2=%.3f \t\\pm %.3f eV \n E0=%.3f \t\\pm %.3f eV \n B=%.3f \t\\pm %.3f mc\n N=%.3f \t\\pm %.3f',...
    chi2min,dof,((M.mnuSq_i+par(1))),err(1),(par(2)),err(2),...
    (M.BKG_RateSec_i+par(3))*1e3,err(3)*1e3,par(4),err(4));

legend([hdata hfit1],'Data',myfitleg,'Location','NorthEast') ;
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Counts/qU','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
title(sprintf('Run %u Data - Samak Fit (%s) \n(18575 eV - qU_{min}) = %.0f eV',RunNr,chi2,M.Q_i-Data(exclDataStart,1)));
PrettyFigureFormat;

subplot(2,1,2);
%  hdata1 = plot(Data(:,1),(Data(:,2)-M.TBDIS)./Data(:,3),...
hdata1 = errorbar(Data(:,1),(Data(:,2)-M.TBDIS),Data(:,3),...
      'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
  
hold on
hdata2 = boundedline(Data(:,1),zeros(26,1),sqrt(diag(FitCM)),'alpha');
%plot(M.qU,zeros(M.nqU,1),'--k');
%plot(M.qU, (Data(:,2)-M.TBDIS)./Data(:,3),'x');
hold off;
xlabel('qU (eV)');
ylabel('Residuals');
grid on;
xlim([min(M.qU(exclDataStart:end)) max(M.qU)]);
PrettyFigureFormat;

if strcmp(saveplot,'ON')
    if  exclDataStart==1 &&exclDataStop == 9999  
        excl = sprintf('all');
    else
        excl = sprintf('excl%u',exclDataStart-1);
    end
    switch rhodScan
        case 'OFF'
    save_name = sprintf('./plots/Run%u/FitRun%u_chi2-%s-%s.eps',RunNr,RunNr,chi2,excl);
        case 'ON'
    save_name = sprintf('./plots/RhoD_Scan/FitRun%u_chi2-%s-%s-%g-molPercm2.eps',RunNr,chi2,excl,WGTS_CD_MolPerCm2_local(rd));     
    end
    publish_figure(fig1,save_name);
    close all;
end
end % end RhoD-Scan
save('./plots/RhoD_Scan/RhoDScan_Sys.mat');
%% rhod plot
switch rhodScan
    case 'ON'
        rhomin = min(WGTS_CD_MolPerCm2_local./D.WGTS_CD_MolPerCm2*100);
        rhomax = max(WGTS_CD_MolPerCm2_local./D.WGTS_CD_MolPerCm2*100);
        
        fig12 = figure(12345); %Endpoint
        set(fig12, 'Units', 'normalized', 'Position', [0.1, 0.1, 1.7 ,0.88]);
        subplot(2,4,1);
        errorbar(WGTS_CD_MolPerCm2_local./D.WGTS_CD_MolPerCm2*100, parRD(2,:),errRD(2,:),...
           's-', 'Color',orange,'LineWidth',1.5);
        xlabel('Column Density (%)');
        ylabel('E_0 - 18575 (eV)');
        xlim([rhomin rhomax]);
        xticks([20:40:140])
        grid on;
        PrettyFigureFormat;             
  
        subplot(2,4,2); % Neutrino mass
        errorbar(WGTS_CD_MolPerCm2_local./D.WGTS_CD_MolPerCm2*100, (M.mnuSq_i+parRD(1,:)),errRD(1,:),...
             's-', 'Color',orange,'LineWidth',1.5);
        xlabel('Column Density (%)')
        ylabel('m^2 (eV)')
        xlim([rhomin rhomax]);
        xticks([20:40:140])
        grid on;        
        PrettyFigureFormat;
        
        subplot(2,4,5); % Background
        errorbar(WGTS_CD_MolPerCm2_local./D.WGTS_CD_MolPerCm2*100, (M.BKG_RateSec_i+parRD(3,:))*1e3,errRD(3,:)*1e3,...
            's-', 'Color',orange,'LineWidth',1.5);
        xlabel('Column Density (%)');
        ylabel(sprintf('BKG (m counts)'));
        xlim([rhomin rhomax]);
        xticks([20:40:140])
        grid on;
        PrettyFigureFormat;
     
        subplot(2,4,6); % Normalization
        errorbar(WGTS_CD_MolPerCm2_local./D.WGTS_CD_MolPerCm2*100, parRD(4,:),errRD(4,:),...
            's-', 'Color',orange,'LineWidth',1.5);
        xlabel('Column Density (%)');
        ylabel('N');
        xlim([rhomin rhomax]);
        xticks([20:40:140])
        grid on;
        PrettyFigureFormat;
        
        subplot(2,4,[3,4,7,8]); % chi2
        plot(WGTS_CD_MolPerCm2_local./D.WGTS_CD_MolPerCm2*100, chi2RD,...
            's-', 'Color',orange,'LineWidth',2);
        xlabel('Column Density (%)');
        ylabel('\chi^2/22 dof');
        xlim([rhomin rhomax]);
        xticks([20:20:150])
        grid on;
        PrettyFigureFormat;
        legend(sprintf('100%% Column Density = %g mol/cm2',D.WGTS_CD_MolPerCm2));
        title(sprintf('Very First Tritium Run %u\nSamak Fit (%s) Results: %.0f eV below E0',RunNr,chi2,M.Q_i-Data(exclDataStart,1)));
       
        if strcmp(saveplot,'ON')
        save_name = sprintf('./plots/RhoD_Scan/RhoDScan-Run%u.eps',RunNr);
        %publish_figure(fig12, save_name);  
        export_fig(save_name)
        end
end
end
end