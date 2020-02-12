%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intro info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../../../Samak2.0'));

% CD = column density

% Fit options
fixPar = '1 11 12';
fitter = 'matlab';
chi2name = 'chi2Stat';
excludeDataStart = 9;
%MRA = [];

% Check if there is already a set of models to skip directly to the fit
if ~exist('MRA','var')


% Arrange runs (datasets) in cells depending on CD
RCD{1} = [40794:40804]; %24
RCD{2} = [40763,40764,40765,40766]; % 48
RCD{3} = [40926:40935]; %71
RCD{4} = [40668]; % 100

% Allocate memory
MRA = cell(4,1);
CD = zeros(4,1);
TimeSecFrac = zeros(26,4);
StackFileName = MRA;
Stacked = MRA;
CovMat = MRA;
MultiModel = MRA;
TBDIS = TimeSecFrac;
qU = TBDIS;

% Loop to extract data
for r = 1:4
    MRA{r} = MultiRunAnalysis('RunList',RCD{r},'chi2',chi2name,...
        'exclDataStart',excludeDataStart,...
        'fixPar','1 5 6','ringCutFlag','ex2','AnaFlag','StackPixel',...
        'DataEffCorr','OFF');
    Stacked{r} = MRA{r}.StackedRuns;
    StackFileName{r} = MRA{r}.StackFileName;
    CD(r) = MRA{r}.RunData.WGTS_CD_MolPerCm2;
    TBDIS(:,r) = MRA{r}.RunData.TBDIS; % counts
    qU(:,r) = MRA{r}.RunData.qU; % retarding potential
    MultiModel{r} = MRA{r}.ModelObj; %one model for each dataset
    CovMat{r} = MRA{r}.FitCM(excludeDataStart:end,excludeDataStart:end);
end

% Build (mega) covariance matrix to use with four datasets
MultiFitCM = blkdiag(CovMat{:});

% Data of all datasets to fit
Data = [qU, TBDIS, sqrt(TBDIS)];

% CD Run Summary = 4.45983560e+17
CDRatio = CD/CD(4);


% First
% Build special class for the fitter class
% computeTBDDS and TBDIS now iterate over the four models 
rhoDM = RHOD('TBDCell',MultiModel);
rhoDM.ComputeTBDDS('N_bias',[0 0 0 0],'B_bias',[0 0 0 0]);
rhoDM.ComputeTBDIS();
end

%function RhoDScan(obj,CD,varargin)
%% Column Density Scan

% p=inputParser;
% p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
% p.parse(varargin{:});
% saveplot = p.Results.saveplot;

saveplot = 'OFF'; % parser parameter

WGTS_CD_MolPerCm2_local = CD(4).*(0.5:0.1:1.5);
%WGTS_CD_MolPerCm2_local = CD(4).*(1.25);

parRD                   = zeros(12,numel(WGTS_CD_MolPerCm2_local));
errRD                   = zeros(12,numel(WGTS_CD_MolPerCm2_local));
chi2RD                  = zeros(numel(WGTS_CD_MolPerCm2_local),1);

for rd = 1:numel(WGTS_CD_MolPerCm2_local)
    
    % Loop on models to rebuild them with the different CDs
    for r = 1:4
        MultiModel{r} = ref_RunAnalysis(StackFileName{r},'ex2','','ISCS','Theory','recomputeRF','OFF',...
            'WGTS_CD_MolPerCm2',CDRatio(r)*WGTS_CD_MolPerCm2_local(rd));
    end
    
    % Second
    % reinitilize new class with the new CD
    rhoDM = RHOD('TBDCell',MultiModel);
    rhoDM.ComputeTBDDS('N_bias',[0 0 0 0],'B_bias',[0 0 0 0]);
    rhoDM.ComputeTBDIS();
    
    
    % Fit
    F = FITC('SO',rhoDM,'DATA',Data,'fitter',fitter,...
        'chi2name',chi2name,...
        'fixPar',fixPar,...
        'exclDataStart',excludeDataStart,...
        'COVMAT',MultiFitCM);
    
    
    FitResult = struct(...
        'par',F.RESULTS{1},....
        'err',F.RESULTS{2},....
        'chi2min',F.RESULTS{3},...
        'errmat',F.RESULTS{4},...
        'dof',F.RESULTS{5});
    %F.catss; F.CATSplot;
    
    % Fit and plot
    if strcmp(saveplot,'ON')
        for mm = 1:4
            FitResultperSpectrum = struct(...
                'par',FitResult.par([1,2,2+mm,6+mm]),....
                'err',FitResult.err([1,2,2+mm,6+mm]),....
                'chi2min',F.RESULTS{3},...
                'errmat',F.RESULTS{4},...
                'dof',F.RESULTS{5});
            P = PLOTC('Xdata',qU(:,mm),'Ydata',TBDIS(:,mm),...
                'ModelObj',MultiModel{mm},'saveplot','export',...
                'titleFlag','RhoDScan','RunList',RCD{mm},'startqU',excludeDataStart,...
                'FitResult',FitResultperSpectrum);
            P.savename = ['RhoD',num2str(mm),'.pdf'];
            names{mm} = ['plots/',P.savename];
            P.PlotSpectrumAndResiduals();
            close all;
        end
        P.savename = 'rhoDFits';
        P.appendPDFex(names);
        
    end
    

    parRD(:,rd) =  FitResult.par;
    errRD(:,rd) =  FitResult.err;
    chi2RD(rd)  =  FitResult.chi2min;
    
end

dofRD       =  FitResult.dof;

%% Find minimum and asymmetric errors
[~,minIndex] = min(chi2RD);
rhoDinPercentage = WGTS_CD_MolPerCm2_local./CD(4)*100;
dataPointsforPlot = [max(1,minIndex-2):min(length(chi2RD),minIndex+2)];
[polyCoeff,structForErrors] = polyfit(rhoDinPercentage(dataPointsforPlot)',chi2RD(dataPointsforPlot),2);

CDminPercentage = -0.5*polyCoeff(2)/polyCoeff(1);
CDmin = CDminPercentage*CD(4)/100;
chi2RDmin = polyval(polyCoeff,CDminPercentage);

chi2err = chi2RDmin + 1;
chi2RDLeftofMin = chi2RD(rhoDinPercentage < CDminPercentage);
rhoDLeftofMin = rhoDinPercentage(rhoDinPercentage < CDminPercentage);
chi2RDRightofMin = chi2RD(rhoDinPercentage > CDminPercentage);
rhoDRightofMin = rhoDinPercentage(rhoDinPercentage > CDminPercentage);
try
    rhoDlowerUnc = interp1(chi2RDLeftofMin,rhoDLeftofMin,chi2err,'spline')*CD(4)/100;
    rhoDupperUnc = interp1(chi2RDRightofMin,rhoDRightofMin,chi2err,'spline')*CD(4)/100;
catch
    try
        rhoDupperUnc = interp1(chi2RDRightofMin,rhoDRightofMin,chi2err,'spline')*CD(4)/100;
        rhoDlowerUnc = 2*CDmin - rhoDupperUnc;
    catch
        rhoDlowerUnc = interp1(chi2RDLeftofMin,rhoDLeftofMin,chi2err,'spline')*CD(4)/100;
        rhoDupperUnc = 2*CDmin - rhoDlowerUnc;
    end
end


%% Plot Results
rhomin = min(WGTS_CD_MolPerCm2_local./CD(4)*100);
rhomax = max(WGTS_CD_MolPerCm2_local./CD(4)*100);

fig12 = figure(12345); %Endpoint
screensize = (get(0, 'Screensize'));
set(fig12, 'Position', [screensize(1:2)*0.95,screensize(3)*0.95,screensize(4)*0.95]);
%set(fig12, 'Units', 'normalized', 'Position', [0.9, 0.9, 1.4, 1.5]);
subplot(2,4,1);
errorbar(WGTS_CD_MolPerCm2_local./CD(4)*100, parRD(2,:),errRD(2,:),...
    's-', 'Color',rgb('CadetBlue'),'LineWidth',1.5);
xlabel('Column Density (%)');
ylabel('E_0 - 18575 (eV)');
xlim([rhomin rhomax]);
xticks([20:40:180])
grid on;
PrettyFigureFormat;

subplot(2,4,2); % Neutrino mass
errorbar(WGTS_CD_MolPerCm2_local./CD(4)*100, (0+parRD(1,:)),errRD(1,:),... % 0 is the initial neutrino mass squared for fit
    's-', 'Color',rgb('CadetBlue'),'LineWidth',1.5);
xlabel('Column Density (%)')
ylabel('m^2 (eV)')
xlim([rhomin rhomax]);
xticks([20:40:180])
grid on;
PrettyFigureFormat;

subplot(2,4,5); % Background
errorbar(WGTS_CD_MolPerCm2_local./CD(4)*100, sum(rhoDM.BKG_RateSec_i+parRD(3:6,:))*1e3/4,sqrt(sum(errRD(3:6,:).^2))*1e3/4,...
    's-', 'Color',rgb('CadetBlue'),'LineWidth',1.5);
xlabel('Column Density (%)');
ylabel(sprintf('BKG (mcps) (average)'));
xlim([rhomin rhomax]);
xticks([20:40:180])
grid on;
PrettyFigureFormat;

subplot(2,4,6); % Normalization
errorbar(WGTS_CD_MolPerCm2_local./CD(4)*100, 1+mean(parRD(7:10,:)),sqrt(sum(errRD(7:10,:).^2))/4,...
    's-', 'Color',rgb('CadetBlue'),'LineWidth',1.5);
xlabel('Column Density (%)');
ylabel('N (average)');
xlim([rhomin rhomax]);
xticks([20:40:180])
grid on;
PrettyFigureFormat;

subplot(2,4,[3,4]) % chi2
hold on
plotchi2 = plot(WGTS_CD_MolPerCm2_local./CD(4)*100, chi2RD,...
    's-', 'Color',rgb('CadetBlue'),'LineWidth',2);
plotchi2fit = plot(WGTS_CD_MolPerCm2_local(dataPointsforPlot)./CD(4)*100,...
    polyval(polyCoeff,WGTS_CD_MolPerCm2_local(dataPointsforPlot)./CD(4)*100),'o-','color','red');
plotchi2fit2 = plot(WGTS_CD_MolPerCm2_local(dataPointsforPlot)./CD(4)*100,...
    polyval(polyCoeff,WGTS_CD_MolPerCm2_local(dataPointsforPlot)./CD(4)*100),'o-','color','red');
hold off
xlabel('Column Density (%)');
ylabel(['\chi^2/',num2str(dofRD),' dof']);
xlim([rhomin rhomax]);
xticks([20:20:180])
set(gca,'YAxisLocation','right');
grid on;
PrettyFigureFormat;
title(sprintf('First Tritium Runs from %u to %u \nSamak Fit (%s) Results: %.0f eV below E0',min(RCD{1}),max(RCD{4}),chi2name,rhoDM.Q_i-rhoDM.TBDCell{1}.qU( excludeDataStart,:)));

legend([plotchi2,plotchi2fit,plotchi2fit2],sprintf('100%% Column Density = %.3g mol/cm2', CD(4)),...
    sprintf('\\rho d_{min} = %0.4g (+%0.4g -%0.4g )', CDmin, rhoDupperUnc- CDmin,-1*( rhoDlowerUnc- CDmin)),...
    sprintf('\\rho d_{min} = %0.4g %% (+%0.4g %% -%0.4g %%)', CDmin/CD(4)*100, (rhoDupperUnc-CDmin)/CD(4)*100,(-1*(rhoDlowerUnc-CDmin))/CD(4)*100),...
    'interpreter','latex');


subplot(2,4,[7,8]) % p-values
plot(WGTS_CD_MolPerCm2_local./CD(4)*100, 1-chi2cdf(chi2RD,FitResult.dof),...
    's-', 'Color',rgb('CadetBlue'),'LineWidth',2);
xlabel('Column Density (%)');
ylabel('p-value');
xlim([rhomin rhomax]);
xticks([20:20:180])
set(gca,'YAxisLocation','right');
grid on;
PrettyFigureFormat;
export_fig(['plots/','rhodScan','ShortRangeWithCovMat'],'-pdf');

if strcmp(saveplot,'ON')
    save_name = sprintf('./plots/RhoD_Scan/RhoDScan%s-%s.pdf',chi2name,'rhoD');
    export_fig(fig12,save_name);
    
end

%end





