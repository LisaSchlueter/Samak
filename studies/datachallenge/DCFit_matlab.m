%
%            FIT INTEGRAL SPECTRUM
%           Fit Active Neutrino Mass


%clear;
tic
addpath(genpath('../../../Samak2.0'));

Nfit = 100; % Number of Fits
mnu_real = 0.0; bck_real = 0.0;
FSD = 'DOSS'; FSDadd = ''; FSDplot = 'only T2';

FitFlag         = 'ON';
savespectrum    = 'OFF';
plotresults     = 'OFF';
plotdist        = 'ON';
saveresults     = 'OFF';
savedist        = 'OFF';
fitter          = 'mat';


% Parametrization: True Value
mnuSq_t = (0.00)^2;

% Loop on fits
npar       = 6;
fit_p      = zeros(Nfit,12);
fit_xi2min = zeros(1,Nfit);
bck_data   = zeros(1,Nfit);

% Init
opt_katrin = {...
    'TD','Flat30',...
    'TimeSec',3*365.25*24*60*60,...
    'Mode','Sim'};

opt_bin = {...
    'nTeBinningFactor',5};

opt_wgts = {...
    'WGTS_CosMaxAAngle',0.6324,...
    'WGTS_Tp',1,...
    'WGTS_DTHTr',1,...
    'WGTS_FTR_cm',4.1,...
    'WGTS_CD_MolPerCm2',5e17,...
    'WGTS_B_T',3.6};

opt_mace = {...
    'MACE_Bmax_T',6,...
    'MACE_Ba_T',9e-4,...
    'MACE_R_eV',2.79};

opt_wgtsmace = {...
    'KTFFlag','SSCW_DCperqU'}; 
                       
opt_fpd = {...
    'FPD_Segmentation','OFF',...
    'FPD_Pixel', 0};

opt_bkg = {...
    'BKG_Flag','ON',...
    'BKG_Type','FLAT',...
    'BKG_RateAllFPDSec',0.0};

opt_fsd= {...
    'TTFSD','DOSS',...
    'DTFSD','OFF',...
    'HTFSD','OFF'};

opt_theocorr = {...
    'ScreeningFlag','OFF',...
    'FiniteExtChargeFlag','OFF',...
    'EEexchangeFlag','OFF',...
    'RecoilCoulombFlag','OFF',...
    'RadiativeFlag','OFF',...
    'RecoilWmVmAFlag','OFF',...
    'WintFiniteSizeFlag','OFF'};

% Tritium spectrum definition
%global A ; 
A = TBD(...
    opt_katrin{:},...
    opt_wgts{:},...
    opt_mace{:},...
    opt_wgtsmace{:},...
    opt_fsd{:},...
    opt_fpd{:},...
    opt_bkg{:},...
    opt_bin{:},...
    opt_theocorr{:},...
    'mnuSq_i',mnuSq_t);

A.ComputeTBDDS();
A.ComputeTBDIS('IStype','SuperFAST');
toc
% progressbar('Neutrino Mass Fits');

% produce and save spectrum
if strcmp(savespectrum,'ON')
    TimePerBin = A.TimeSec*A.qUfrac;
    for ii = 1:1
        A.ComputeTBDIS('IStype','SIMP');
        %A.AddStatFluctTBDIS();
        TBDIS = real([A.qU,A.TBDIS./TimePerBin,A.TBDISE./TimePerBin]);
        save(['DC_Data/dtdc_v1_samak_','nostatfluct','2.txt'],'TBDIS','-ascii');
        %save(['DC_Data/RFperqU/test_samak_','nostatfluct','.mat'],'TBDIS','-v7.3','-nocompression')
        disp(ii)
    end
end


if strcmp(FitFlag,'ON')
    for f = 1:Nfit
        tic
        % Initializing Fit
        DataName = sprintf('DC_Data/FitDC/K3y_flat30_Doss01T2_%04d.txt', f);
                %DataName = sprintf('DC_Data/RFperqU/test_samak%04d.mat', f);
        % DataName = 'DC_Data/dtdc_v1_ssc_nostatfluct.txt';
        DataWithHeaders = importdata(DataName);
        %DataWithHeaders = load(['DC_Data/RFperqU/test_samak_',num2str(f),'.mat']);
        %         [nqU_data,columns_data] = size(DataWithHeaders.data);
        
        %Data = importdata(DataName);
        %         Data = [Data sqrt(Data(:,2).*(A.TimeSec*A.qUfrac))./(A.TimeSec*A.qUfrac)];
        %        Data = load(DataName);
        %       Data = Data.TBDIS;
        %D = load('runmat1.mat');

        %Data = D.TBDISallPixels(:,1:3);
        
        Data = [DataWithHeaders.data(:,1)+18575,...
            DataWithHeaders.data(:,2),...
            DataWithHeaders.data(:,3)];
        
        
        
        % Init
        i_mnu = 0.00^2;
        i_Q = 0;
        i_B = Data(end,2);
        i_N = Data(1,2)/(A.TBDIS(1)/(A.qUfrac(1)*A.TimeSec)+i_B) - 1;
        i_m4 = 0;
        i_s4 = 0;
        
        parnames = ['mSq dQ B dN m4 s4'];
        ParIni = [i_mnu i_Q i_B i_N i_m4 i_s4];
        tmparg = sprintf(['set pri -10;'...
            'fix  5 6 ; min']);
        DataTBD = {Data,A};
        
        switch  fitter
            
            case 'min'
                
                Args = {ParIni, DataTBD, '-c',tmparg};
                [par, err, chi2min, errmat] = fminuit('NuMassChi2DC',Args{:});
                %[par, err, chi2min, errmat] = pminuit('NuMassChi2DC',Args);
            case 'mat'
                options = optimoptions('fminunc','Algorithm','quasi-newton',...
                     'OptimalityTolerance',1e-8,'StepTolerance',1e-9,'FiniteDifferenceType','central');
                TBDfun = @(xx) NuMassChi2DC(xx,DataTBD);
                [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,ParIni,options);
                if exitflag ~= 1
                    [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,par,options);      
                end
                Hessian = 0.5*Hessian(1:4,1:4);
                varcov = inv(Hessian);
                err = sqrt(diag(varcov));
        end
        
        mnuSq_report = A.mnuSq_i+par(1);
        mnu_report = sqrt(abs(mnuSq_report));
        err_mnu = sqrt(abs((mnuSq_report + err(1)))) - mnu_report;
        BCK_fit = A.BKG_RateSec_i+par(3);
        fprintf('===============================================\n');
        fprintf('  m^2       = %g ± %g eV^2\n', mnuSq_report,err(1));
        fprintf('  m         = %g ± %g eV\n', mnu_report,err_mnu);
        fprintf('  dQ        = %g ± %g eV\n', par(2),err(2));
        fprintf('  B         = %g ± %g cps\n', BCK_fit,err(3));
        fprintf('  dN        = %g ± %g\n', par(4),err(4));
        fprintf('  Chi2/dof  = %g/%g\n', chi2min, A.nqU-4);
        fprintf('===============================================\n');
        
        % 1 mnu^2 2 mnu^2_err 3 mnu 4 mnu_err 5 Q 6 Q err
        % 7 bck 8 bck err 9 norm 10 norm err 11 chi2 12 dof
        fit_p(f,:) = [mnuSq_report err(1) mnu_report err_mnu par(2) err(2) ...
            par(3) err(3) par(4) err(4) chi2min A.nqU-4];
                  
        %% Plot Results
        if strcmp(plotresults,'ON')
            FitIS = NuMassModel4parDC(par,A);
            figure(1)
            subplot(2,1,1)
            hold on
            hfit = plot(Data(:,1)-A.Q_i,FitIS,...
                'LineWidth',1,'LineStyle','-','Color','Black');
            hfit2 = plot(Data(:,1)-A.Q_i,FitIS./(A.TimeSec*A.qUfrac),...
                'LineWidth',1,'LineStyle','-','Color','Black');
            hdata = errorbar(Data(:,1)-A.Q_i,Data(:,2),Data(:,3),...
                'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
            hold off
            grid on
            xlabel('qU-E_0 (eV)','FontSize',14);
            ylabel('Counts','FontSize',14);
            title(sprintf('KATRIN - %s - %gy ',A.TD,3));
            set(gca,'FontSize',12);
            mydata = 'Data: m^2 = ??? eV^2';
            myfit2 = sprintf('Fit: m^2 = %.3f \\pm %.3f eV^2\n',mnuSq_report,err(1));
            mychi2 = sprintf('\\chi^2/dof = %.1f/%.0f\n',chi2min,A.nqU-4);
            lgd_fit = legend([hdata hfit hfit2],mydata,myfit2,mychi2,'Location','Best') ;
            lgd_fit.FontSize = 8;
            axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 min(Data(:,2)) max(Data(:,2))*1.1])
            
            subplot(2,1,2)
%             parnomix = par; parnomix(1) = 0; %- A.mnuSq;
%             hfit = plot(Data(:,1)-A.Q_i,...
%                 (Data(:,2))./-1,...
%                 'Color','Black','LineWidth',1,'LineStyle','-');
%             hold on;
%             hdata = errorbar(Data(:,1)-A.Q,...
%                 Data(:,2)./(NuMassModel4par(par,A)./(A.TimeSec*A.qUfrac))-1,...
%                 Data(:,3)./(NuMassModel4par(par,A)./(A.TimeSec*A.qUfrac)),...
%                 'ks','MarkerSize',5,'MarkerFaceColor',0.9*[1 1 1],'Color','Black','LineWidth',1);
            hdata = scatter(Data(:,1)-A.Q,...
                (Data(:,2)-(NuMassModel4par(par,A)./(A.TimeSec*A.qUfrac)))./Data(:,3),...
                'ks','LineWidth',1);
            %             hold off;
            grid on
            xlabel('qU-E_0 (eV)','FontSize',14);
            ylabel('Residuals','FontSize',14);
            set(gca,'FontSize',12);
            %             lgd_res = legend([hdata hfit],'Data/No mass - 1','Fit/No mass - 1',...
            %                 'Location','SouthWest');
            lgd_res = legend([hdata],'(Data-Fit)/uncertainty - 1',...
                'Location','Best');
            lgd_res.FontSize = 8;
            %             axis([min(A.qU-A.Q_i) max(A.qU-A.Q_i)+1 ...
            %                 min((Data(:,2)./(NuMassModel4par(parnomix,A)./(A.TimeSec*A.qUfrac))-1)) ...
            %                 max(Data(:,2)./(NuMassModel4par(parnomix,A)./(A.TimeSec*A.qUfrac)))-1]);
%                         axis([min(A.qU-A.Q_i), max(A.qU-A.Q_i)+1, -13e-4, 13e-4]);
            if strcmp(saveresults,'ON')
                %print(['DC_Data/FigsDC/samakfitSSC',num2str(f)],'-dpng')
                print(['DC_Data/RFTest/samakfitSSC',num2str(f)],'-dpng')
            end
            
            % progressbar(f/Nfit);
        end
        disp(f)
        toc
    end
end
elapsedtime = toc; 

%% Plot distribution
% WORKING ON THIS
if strcmp(plotdist,'ON')
    
    %load('fitDC1000.mat'); 
%     fit_p = zeros(Nfit,11);
%     DATASSC = importdata('step2_optional_by_ssc.txt');
%     fit_p(:,1) = DATASSC.data(:,2); %mnu^2
%     fit_p(:,5) = DATASSC.data(:,4); %E_0 
%     fit_p(:,7) = DATASSC.data(:,6); %bck
%     fit_p(:,9) = DATASSC.data(:,8); %norm
%     fit_p(:,11) = DATASSC.data(:,10); %chi^2
    
    figure(2);
    [mu_mnu,sigma_mnu] = normfit(fit_p(:,1));
    subplot(2,2,1);
    hold on
    [tbins,Nbins,Wbins] = nhist(fit_p(:,1));
    x_mnu = linspace(min(fit_p(:,1)),max(fit_p(:,1)),length(Nbins)*A.nTeBinningFactor);
    binWidth = Wbins(2) - Wbins(1);
    plot(x_mnu,Nfit*binWidth*normpdf(x_mnu,mu_mnu,sigma_mnu));
    hold off
    annotation('textbox',[0.36 0.7 0.2 0.2],'String',...
        {['\sigma = ',num2str(sigma_mnu),' eV^2'],['mean = ',num2str(mu_mnu),' eV^2']},'FitBoxToText','on');
    title('m_\nu^2 [eV^2]');
    xlabel('m_\nu^2 [ev^2]');
    
    [mu_Q,sigma_Q] = normfit(fit_p(:,5));
    subplot(2,2,2);
    hold on
    [tbins,Nbins,Wbins] = nhist(fit_p(:,5));
    x_Q = linspace(min(fit_p(:,5)),max(fit_p(:,5)),length(Nbins)*A.nTeBinningFactor);
    binWidth = Wbins(2) - Wbins(1);
    plot(x_Q,Nfit*binWidth*normpdf(x_Q,mu_Q,sigma_Q));
    hold off
    annotation('textbox',[0.8 0.7 0.2 0.2],'String',...
        {['\sigma = ',num2str(sigma_Q),' eV'],['mean = ',num2str(mu_Q),' eV']},'FitBoxToText','on');
    title('Q [eV]');
    xlabel('Q [eV]');
    
    [mu_b,sigma_b] = normfit(fit_p(:,7) - bck_real);
    subplot(2,2,3);
    hold on
    [tbins,Nbins,Wbins] = nhist(fit_p(:,7) - bck_real);
    x_b = linspace(min(fit_p(:,7) - bck_real),max(fit_p(:,7) - bck_real),length(Nbins)*A.nTeBinningFactor);
    binWidth = Wbins(2) - Wbins(1);
    plot(x_b,Nfit*binWidth*normpdf(x_b,mu_b,sigma_b));
    hold off
    annotation('textbox',[0.36 0.2 0.2 0.2],'String',...
        {['\sigma = ',num2str(sigma_b),' cps'],['mean = ',num2str(mu_b),' cps']},'FitBoxToText','on');
    title('background [cps]');
    xlabel('background [cps]');
    
    [mu_N,sigma_N] = normfit(fit_p(:,9));
    subplot(2,2,4);
    hold on
    [tbins,Nbins,Wbins] = nhist(fit_p(:,9));
    x_N = linspace(min(fit_p(:,9)),max(fit_p(:,9)),length(Nbins)*A.nTeBinningFactor);
    binWidth = Wbins(2) - Wbins(1);
    plot(x_N,Nfit*binWidth*normpdf(x_N,mu_N,sigma_N));
    hold off
    annotation('textbox',[0.8 0.2 0.2 0.2],'String',...
        {['\sigma = ',num2str(sigma_N),''],['mean = ',num2str(mu_N),'']},'FitBoxToText','on');
    title('normalization');
    xlabel('normalization');
    
    set(get(axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off'),'Title'),'Visible','on')
    title(['Fit ',FSD,' ',FSDplot,' FSD']);
    set(gcf,'units','pixels','position',[1,41,1536,748.8])
    
    if strcmp(savedist,'ON')       
        print(['DC_Data/RFTest/DCDist',FSD,FSDadd],'-dpng')
    end
    
    figure(3)
    hold on
    [tbins,Nbins,Wbins] = nhist(fit_p(:,11));
    x_chi2 = linspace(min(fit_p(:,11)),max(fit_p(:,11)),length(Nbins)*A.nTeBinningFactor);
    binWidth = Wbins(2) - Wbins(1);
    plot(x_chi2,Nfit*binWidth*chi2pdf(x_chi2,A.nqU-4))
    hold off
    annotation('textbox',[0.7 0.7 0.2 0.2],'String',['DoF = ',num2str(A.nqU-4),''],'FitBoxToText','on');
    title(['\chi^2 Fit ',FSD,' ',FSDplot,' FSD']);
    xlabel('\chi^2');
    if strcmp(savedist,'ON')
        print(['DC_Data/RFTest/DCChi2',FSD,FSDadd],'-dpng')
    end
    
    CorrPlot = corrplotm(fit_p(:,[1 5 7 9]),'varNames',{'m_\nu^2','Q','BCK','NORM'});
    ax = axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
    set(get(ax,'Title'),'Visible','on');
    ax.Title.String = {' ',[FSD,' ',FSDplot,' FSD']};
    set(gcf,'units','pixels','position',[1,41,1536,748.8])
    if strcmp(savedist,'ON')
        print(['DC_Data/RFTest/CorrPlot',FSD,FSDadd],'-dpng')
    end
%     set(gcf,'units','pixels','position',[1,41,1536,748.8])
end

