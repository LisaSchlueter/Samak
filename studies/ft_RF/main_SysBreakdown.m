%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Study of Response Function Systematics
% Breakdown of ELoss - BFields - Column Density / Cross Section Part
% CovMatrices and Fit
% L.Schlueter 05/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('../../../Samak2.0'));

plotFlag = 'ON';
CMFlag = 'Frac'; %'CM'
TimeSec = 24*60*60*3;
mtd = 'FT-TL4';
Fitter = 'Minuit';
nfits = 500;

%% Build/Load Covariance Matrices
A = InitKatrin_ft_RF('TD', mtd, 'TimeSec', TimeSec);
myEffects = struct(...
    'RF_EL','ON',...  % Response Function(RF) EnergyLoss
    'RF_BF','ON',...  % RF B-Fields
    'RF_RX','ON');  % RF Column Density, Cross Section
C = CovarianceMatrix('StudyObject',A, 'nTrials',1000,'SysEffect',myEffects ,'RecomputeFlag','OFF','SanityPlots','OFF');
C.ComputeCM_RF; %all ON
%C.PlotCM('saveplot',plotFlag,'savename','RF-all');

C.SysEffect.RF_BF = 'ON';
C.SysEffect.RF_RX = 'OFF';
C.SysEffect.RF_EL = 'OFF';
C.ComputeCM_RF;
%C.PlotCM('saveplot',plotFlag,'savename','RF-BF');

C.SysEffect.RF_BF = 'OFF';
C.SysEffect.RF_RX = 'ON';
C.SysEffect.RF_EL = 'OFF';
C.ComputeCM_RF;
%C.PlotCM('saveplot',plotFlag,'savename','RF-RX');

C.SysEffect.RF_BF = 'OFF';
C.SysEffect.RF_RX = 'OFF';
C.SysEffect.RF_EL = 'ON';
C.ComputeCM_RF;
%C.PlotCM('saveplot',plotFlag,'savename','RF-EL');

C.ComputeCM_STAT;
%% Fit
A = InitKatrin_ft_RF('TD', mtd, 'TimeSec', TimeSec); % Model
D = InitKatrin_ft_RF('TD', mtd, 'TimeSec', TimeSec); % Data (with Bias)

myCovMatFrac = {zeros(A.nqU,A.nqU), C.MultiCovMatFrac.CM_RF, C.MultiCovMatFrac.CM_RF_BF, C.MultiCovMatFrac.CM_RF_EL, C.MultiCovMatFrac.CM_RF_RX};
myCovMat     = {zeros(A.nqU,A.nqU), C.MultiCovMat.CM_RF, C.MultiCovMat.CM_RF_BF, C.MultiCovMat.CM_RF_EL, C.MultiCovMat.CM_RF_RX};

Fitpar = zeros(6,length(myCovMatFrac),nfits); %Save Results
Fiterr = zeros(6,length(myCovMatFrac),nfits);
FitE0mean = zeros(length(myCovMatFrac),1);
FitE0errmean = zeros(length(myCovMatFrac),1);
Fitchi2min = zeros(length(myCovMatFrac),nfits);
FitLabel = {'Stat', 'RF','RF-BF', 'RF-EL','RF-RX'};
tic
for i=1:length(myCovMatFrac)
        for f=1:nfits %fits with stat fluct
            D.ComputeTBDDS(); D.ComputeTBDIS('IStype','SIMP');
            D.AddStatFluctTBDIS; % Fluctuate the Data: stat
            Data = [D.qU,D.TBDIS,sqrt(D.TBDIS)];
            FitCovMatFrac = C.MultiCovMatFrac.CM_STAT + myCovMatFrac{i};
            FitCovMat = C.MultiCovMat.CM_STAT + myCovMat{i};
            
            switch CMFlag
                case 'Frac'
                    DataTBD = {Data,A, FitCovMatFrac};  %Data, Modell, Covariance Matrix
                case 'CM'
                    DataTBD = {Data,A, FitCovMat};  %Data, Modell, Covariance Matrix
            end
            % Init Fit Parameters
            i_mnu      = 0; i_Q        = 0;
            i_B        = D.TBDIS(end)./(D.TimeSec*D.qUfrac(end));
            i_N        = 0;
            i_m4       = 0; i_s4       = 0;          
            ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4];
            
            switch Fitter
                case 'Minuit'
                    parnames = ['mSq dQ B dN m4 s4'];
                    tmparg = sprintf(['set pri -1 ; fix  5 6 ; min ; ' ]);
                    Args = {ParIni, DataTBD, '-c',tmparg};
                    switch CMFlag
                        case 'Frac'
                            [par, err, chi2min, errmat] = fminuit('NuMassChi2E0_CovFrac',Args{:});
                        case 'CM'
                            [par, err, chi2min, errmat] = fminuit('NuMassChi2E0_Cov',Args{:});
                    end
                case 'Matlab'                   
                options = optimoptions('fminunc','Algorithm','quasi-newton',...
                    'OptimalityTolerance',1e-7,'StepTolerance',1e-7,'FiniteDifferenceType','central');
                TBDfun = @(xx) NuMassChi2E0_Cov(xx,DataTBD);
                [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,ParIni,options);
                if exitflag ~= 1
                    [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,par,options);
                end
                errmat = 0.5*Hessian;
                varcov = inv(errmat);
                err = sqrt(diag(varcov));
            end
            %save results
            Fitpar(:,i,f) = par;
            Fiterr(:,i,f) = err;
            Fitchi2min(i,f) = chi2min;
           
            fprintf(2,'----------------------------------------------\n');
            fprintf('  Processing \t= %g %% \n',f/nfits*100);
            mnuSq_report = A.mnuSq_i+par(1);
            mnu_report = sqrt(abs(mnuSq_report));
            err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
            BCK_fit = A.BKG_RateSec_i+par(3);
            fprintf('===============================================\n');
            fprintf('  m^2       = %g +- %g eV^2\n', mnuSq_report,err(1));
            fprintf('  m         = %g +- %g eV\n', mnu_report,err_mnu);
            fprintf('  dQ        = %g +- %g eV\n', par(2),err(2));
            fprintf('  B         = %g +- %g cps\n', BCK_fit,err(3));
            fprintf('  dN        = %g +- %g\n', par(4),err(4));
            fprintf('  Chi2/DoF  = %g/%g\n', chi2min, A.nqU-4);
            fprintf('===============================================\n');            
        end 
        FitE0mean = mean(Fitpar(2,i,:));
        FitE0errmean = mean(Fiterr(2,i,:));
        %%
        f1 = figure(1);
        nhist(Fitchi2min(i,:));
        xlabel('\chi^2/14 dof');
        ylabel('');
        legend(sprintf('nfits = %u \nmean =%.2f/ %u dof \nsigma=%.2f',nfits, mean(Fitchi2min(i,:)),A.nqU-4,sqrt(var(Fitchi2min(i,:)))));     
        fig_name = sprintf('./plots/SysBreakdown/Fit_%s-%ufits_chi2.pdf',FitLabel{i}, nfits);
        publish_figure(f1,fig_name);
        
        f11= figure(2);
        nhist(Fitpar(2,i,:));
        xlabel('(E0)_{eff} - 18575 [eV]')
        ylabel('');
        legend(sprintf('nfits =%u \nmean =%.2f eV\nsigma=%.2f eV',nfits, mean(Fitpar(2,i,:)),sqrt(var(Fitpar(2,i,:)))));
        fig_name = sprintf('./plots/SysBreakdown/Fit_%s-%ufits_E0.pdf',FitLabel{i},nfits);
        publish_figure(f11,fig_name);
        
        
end
toc
x = linspace(1,length(myCovMatFrac),length(myCovMatFrac));
f3 = figure(3);
subplot(2,1,1);
bar(x,FitE0mean);
xticklabels({'Stat.'; 'RF'; 'R-BF';'RF-EL'; 'RF-RX'})
ylabel('<E0_{eff} - 18575> [eV]');
title(sprintf('First Tritium Nominal Settings - Endpoint Fit \nMTD %s - %.0fh Run - Response Function CovMat',mtd,TimeSec/(60*60)));
subplot(2,1,2);
bar(x,FitE0errmean);
xticklabels({'Stat.'; 'RF'; 'R-BF';'RF-EL'; 'RF-RX'})
ylabel('<\Delta E0_{eff}> [eV]');
PrettyFigureFormat;
fig_name = sprintf('./plots/SysBreakdown/SysBreakdownFit_-%ufits.pdf',nfits);
publish_figure(f3,fig_name);
%% plot 1 : STACKPLOT
if strcmp(plotFlag,'ON')
f1 = figure(999); %stack plot
C.StudyObject.ComputeTBDDS; C.StudyObject.ComputeTBDIS;
BFieldErr = sqrt(diag(C.MultiCovMat.CM_RF_BF))./C.StudyObject.TBDIS;
ELossErr = sqrt(diag(C.MultiCovMat.CM_RF_EL))./C.StudyObject.TBDIS;
RhoXErr =  sqrt(diag(C.MultiCovMat.CM_RF_RX))./C.StudyObject.TBDIS;
AllErrsum = sqrt(diag(C.MultiCovMat.CM_RF_BF)+diag(C.MultiCovMat.CM_RF_EL)+diag(C.MultiCovMat.CM_RF_RX))./C.StudyObject.TBDIS;
Allin1 =  sqrt(diag(C.MultiCovMat.CM_RF))./C.StudyObject.TBDIS;

b1 = bar(C.StudyObject.qU*1e-03, [ RhoXErr  BFieldErr  ELossErr ],'stack');
mycolors = {rgb('Teal'); rgb('DarkOliveGreen'); rgb('YellowGreen')};
% mycolors = {rgb('Maroon'); rgb('Tomato'); rgb('Salmon')};
set(b1,{'FaceColor'},mycolors);
hold on;
%bar(C.StudyObject.qU*1e-03, C.StudyObject.TBDISE./C.StudyObject.TBDIS, 'FaceAlpha', 0);

% bar(C.StudyObject.qU*1e-03, AllErrsum, 'FaceAlpha',0);
bar(C.StudyObject.qU*1e-03, Allin1, 'FaceAlpha',0);

legend('RhoXsection: 10% - 1%' ,'BField: 1%','ELoss 1%', 'Stat. (24 hours)', 'Location','northwest');%'Combined','All in 1',
xlabel('qU (keV)', 'FontSize',15);
ylabel('Rel. Uncertainty on Integral Spectrm', 'FontSize',15);
set(gca, 'XScale', 'log');
PrettyFigureFormat;
hold off;
title(sprintf('First Tritium Nominal Settings \n- Response Function Uncorrelated Uncertainties %u Trials-', C.nTrials));

s_file = sprintf('./plots/SysBreakdown/SysBreakdownRF_%s',A.TD);
%publish_figure(f1,strcat(s_file,'.eps'));
%export_fig(strcat(s_file,'.fig'));
            
%% plot 2: Covariance Matrices

f2 = figure(222);
set(f2, 'Units', 'normalized', 'Position', [0.1, 0.1, 1 ,1]);
subplot(2,2,1);
imagesc(C.MultiCovMatFrac.CM_RF);
colorbar;
colormap(flipud(gray));
title('Fractional CM','FontSize', 10);
PrettyFigureFormat;
pbaspect([1 1 1])

subplot(2,2,2);
imagesc(C.MultiCovMatFrac.CM_RF_RX);
colorbar;
title('RhoXSection(10% - 1%)','FontSize', 10);
PrettyFigureFormat;
pbaspect([1 1 1])

subplot(2,2,3);
imagesc(C.MultiCovMatFrac.CM_RF_EL);
colorbar;
title('ELoss (1%)','FontSize', 10);
PrettyFigureFormat;
pbaspect([1 1 1])

subplot(2,2,4);
imagesc(C.MultiCovMatFrac.CM_RF_BF);
colorbar;
title('BField (1%)','FontSize', 10);
PrettyFigureFormat;
pbaspect([1 1 1])

cm_file = sprintf('./plots/SysBreakdown/CM_RF_%s_%u',A.TD, C.nTrials);
publish_figure(f2,strcat(cm_file,'.eps'));
export_fig(strcat(cm_file,'.fig'));
end



