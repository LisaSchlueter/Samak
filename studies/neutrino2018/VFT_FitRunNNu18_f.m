function [parA , err , chi2min , dof] =VFT_FitRunNNu18_f(varargin)

addpath(genpath('../../../Samak2.0'));


    % Parser
    p = inputParser;
    p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('RunNr',40257,@(x)isfloat(x) && x>=0);
    p.addParameter('exclDataStart',7,@(x)isfloat(x) && x>=0);
    p.addParameter('chi2','CM',@(x)ismember(x,{'CM','Stat'}));
    p.parse(varargin{:});
    display           =    p.Results.display;
    RunNr             =    p.Results.RunNr;
    exclDataStart     =    p.Results.exclDataStart; % Start Fit at bin i (i=1 -> qU(min)=-1.7keV, i=9 ->qU(min)=200eV, i=7 ->qU(min)=400eV)
    chi2              =    p.Results.chi2;

%% Read Data
D = importdata([num2str(RunNr),'.mat']);
switch RunNr
    case {40259,40260,40263, 40264, 40265,40266}
        ErrCount = sqrt(D.TBDIS(1:end));
        if RunNr == 40260
            ErrCount(end)=999;
        end
        Data = [D.qU, D.TBDIS, ErrCount];
    case {40257,40258}
        Data = [D.qU(3:end), D.TBDIS(3:end), sqrt(D.TBDIS(3:end))];
end
%% Model
M = ref_runsummaries(RunNr,'ISCS','Theory','recomputeRF','OFF');
M.ComputeTBDDS; M.ComputeTBDIS;

%% Covariance Matrix
cm = importdata('CovMat_Run40257_0.1rhod_0.02BField_0.02ISX_0.03FSD_0.023TASR_TC.mat');
covmatfrac = cm.CovMatFracCombi;
CM_Stat = diag(Data(:,3).^2);
CMFrac_Stat = diag(1./Data(:,2));
FitCM = Data(:,2).*covmatfrac.*Data(:,2)' + CM_Stat;
FitCMFrac = covmatfrac + CMFrac_Stat; 

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

par  = F.RESULTS{1}; 
parA = [par(1) M.Q+par(2) M.BKG_RateSec_i+par(3) par(4)]
err = F.RESULTS{2};
chi2min = F.RESULTS{3};
dof = numel(M.qU(exclDataStart:end))-4;

%% Plot Spectrum + Fit
switch display
    case 'ON'
fig1 = figure(5);
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,1]);
subplot(2,1,1);
hfit1 = boundedline(M.qU,M.TBDIS,diag(sqrt(FitCM)),'alpha');
hold on;
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold off;
set(gca, 'YScale', 'log')
xlabel('qU (eV)');
ylabel('counts')
xlim([min(M.qU(exclDataStart:end)) max(M.qU)]);
myfitleg = sprintf(' Fit: \\chi2 / dof=%.2f/%.0f \n m^2=%.2f \t\\pm %.2f eV \n E0=%.2f \t\\pm %.2f eV \n B=%.2f \t\\pm %.2f mcps\n N=%.2f \t\\pm %.2f',...
    chi2min,dof,((M.mnuSq_i+par(1))),err(1),(par(2)),err(2),...
    (M.BKG_RateSec_i+par(3))*1e3,err(3)*1e3,par(4)+1,err(4));
legend([hdata hfit1],'Data',myfitleg,'Location','NorthEast') ;
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Counts/qU','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
title(sprintf('Run %u Data - Samak Fit (%s) \n(18575 eV - qU_{min}) = %.0f eV',RunNr,chi2,M.Q_i-Data(exclDataStart,1)));
PrettyFigureFormat;

subplot(2,1,2);
hdata1 = errorbar(Data(:,1),(Data(:,2)-M.TBDIS),Data(:,3),...
      'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
hdata2 = boundedline(Data(:,1),zeros(26,1),sqrt(diag(FitCM)),'alpha');
hold off;
xlabel('qU (eV)');
ylabel('Residuals');
grid on;
xlim([min(M.qU(exclDataStart:end)) max(M.qU)]);
PrettyFigureFormat;
end

