function [parA, err, chi2min, ndof ] = fit_tbdpixel(varargin)
%
% Fit TBD single Pixel Simulation
%
% Th. Lasserre - CEA Saclay
% Last Updated: Feb. 12 2018
%
%add path
%%cd /home/lisa/Dokumente/Studium/Masterstudium/Masterarbeit/MySamak/samak/studies/tbdallpixels;
addpath(genpath('../../../Samak2.0'));

% Initialization% path

clear par ; clear parnomix;

% Parser
p = inputParser;
p.addParameter('Run',1,@(x)isfloat(x));
p.addParameter('Fitter','Minuit',@(x)ismember(x,{'Minuit','Matlab'}));
p.addParameter('Chi2Type','G',@(x)ismember(x,{'G','P'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('displayPlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('PixelID',1,@(x)isfloat(x) && x>0);
p.addParameter('mnuSq_t',(0.)^2,@(x)isfloat(x)); 
p.addParameter('StatFluct','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('MyMode','Simulation',@(x)ismember(x,{'Bootcamp','Simulation'}));
p.parse(varargin{:});
Run               =    p.Results.Run;
Fitter            =    p.Results.Fitter;
Chi2Type          =    p.Results.Chi2Type;
display           =    p.Results.display;
displayPlot       =    p.Results.displayPlot;
PixelID           =    p.Results.PixelID;
mnuSq_t           =    p.Results.mnuSq_t;
StatFluct         =    p.Results.StatFluct;
MyMode              =    p.Results.MyMode;
%save = ('myfile.mat');
switch MyMode
    case 'Bootcamp'
        switch Run
            case 1
                S = init_bootcampR1('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat1.mat');
                figprefix = 'bc_fit_run1';
            case 2
                S = init_bootcampR2('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat2.mat');
                figprefix = 'bc_fit_run2';
            case 3
                S = init_bootcampR3('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat3.mat');
                figprefix = 'bc_fit_run3';
            case 4
                S = init_bootcampR4('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat4.mat');
                figprefix = 'bc_fit_run4';
            case 5
                S = init_bootcampR5('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat5.mat');
                figprefix = 'bc_fit_run5';
            case 6
                S = init_bootcampR6('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat6.mat');
                figprefix = 'bc_fit_run6';
            case 7
                S = init_bootcampR7('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat7.mat');
                figprefix = 'bc_fit_run7';
            case 8
                S = init_bootcampR8('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat8.mat');
                figprefix = 'bc_fit_run8';
            case 9
                S = init_bootcampR9('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat9.mat');
                figprefix = 'bc_fit_run9';
            case 10
                S = init_bootcampR10('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat10.mat');
                figprefix = 'bc_fit_run10';
            case 110
                S = init_bootcampR110('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat1_10.mat');
                figprefix = 'bc_fit_run1_10';
            case 35410
                S = init_bootcampR35410('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat35410.mat');
                figprefix = 'bc_fit_run35410';
            case 35411
                S = init_bootcampR35411('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat35411.mat');
                figprefix = 'bc_fit_run35411';
            case 35412
                S = init_bootcampR35412('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat35412.mat');
                figprefix = 'bc_fit_run35412';
            case 35413
                S = init_bootcampR35413('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat35413.mat');
                figprefix = 'bc_fit_run35413';
            case 35420
                S = init_bootcampR35420('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat35420.mat');
                figprefix = 'bc_fit_run35420';
            case 35422
                S = init_bootcampR35422('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat35422.mat');
                figprefix = 'bc_fit_run35422';
            case 3541035422
                S = init_bootcampR3541035422('PixelID',PixelID,'mnuSq_i',mnuSq_t);
                D = load('data/runmat3541035422.mat');
                figprefix = 'bc_fit_run354103_5422';
        end
        S.TimeSec = D.TimeSec;
    case 'Simulation'
        figprefix = 'bc_fit_runsim';
        S = init_tbdallpixels('PixelID',PixelID,'mnuSq_i',mnuSq_t);
end
S.mnuSq_i=mnuSq_t; 
%SetLineInitParallPixels(A,'Random','OFF')
%S.BKG_RateSec_i = S.BKG_RateSecallPixels(PixelID);

% Loop on fits
npar       = 6; parA = [0 0 0 0 0 0];
ndof       = S.nqU - 4;

switch Fitter
    case 'Minuit'
    case 'Matlab'
        tbdisf = @(mb,e0b,bb,nb,qu) S.ComputeTBDISf(qu,mb,e0b,bb,nb,0,0);
        tmp = @(mb,e0b,bb,nb,e) interp1(S.qU,tbdisf(mb,e0b,bb,nb,S.qU),e);
        myf = fittype(@(mb,e0b,bb,nb,qu) tmp(mb,e0b,bb,nb,qu),...
            'independent','qu','coefficients', {'mb','e0b','bb','nb'});
end
 
switch MyMode
    case 'Bootcamp'
        quIndex      = (PixelID-1)*3+1;  %[1:3:148*3];
        countIndex   = (PixelID-1)*3+2; %[2:3:148*3];
        countErrIndex = (PixelID-1)*3+3; %[3:3:148*3];
        S.qU     = D.TBDISallPixels(:,quIndex);
        S.qUfrac = (D.qUfrac);
        S.SetTBDDSBinning();
        Counts      = D.TBDISallPixels(:,countIndex);
        CountErrors = D.TBDISallPixels(:,countErrIndex); CountErrors(CountErrors==0)=1;
        Data = [...
            S.qU,...
            Counts,...
            CountErrors];
        S.ComputeNormFactorTBDDS;
        S.BKG_RateSec_i = Data(end,2)./S.TimeSec./S.nqU;
        S.ComputeTBDDS; S.ComputeTBDIS; 
        %S.PlotTBDIS;
        %plot(S.qU,Data(:,2));return;
    case 'Simulation'
        S.BKG_RateSec_i = S.BKG_RateAllFPDSec./S.nPixels;
        S.ComputeTBDDS(); S.ComputeTBDIS(); 
        switch StatFluct
            case 'ON'
                AddStatFluctTBDIS(S);
        end
        Data = [S.qU,S.TBDIS,S.TBDISE];
end


% Random Init
i_mnuSq    = 0;%(0.2*randn)^2;
i_E0       = 0;%0.05*randn;
i_B        = 0;%0.1*S.BKG_RateSec_i*randn;
i_N        = 0;%1e-3*randn;

switch Fitter
    case 'Minuit'
        ParIni = [i_mnuSq i_E0 i_B i_N];
        parnames = ['mSq dE0 B dN'];
        tmparg = sprintf([' set pri -1 ; min ']);
        DataL = {Data,S};
        Args = {ParIni, DataL,'-n',parnames, '-c',tmparg};
        switch Chi2Type
            case 'G'
                [par, err, chi2min, errmat] = pminuit('tbd_chi2',Args);
            case 'P'
                [par, err, chi2min, errmat] = pminuit('tbd_chi2Poisson',Args);
        end
        par = [par 0 0]; err = [err 0 0];
        drawnow('update');
    case 'Matlab'
        opts = fitoptions('StartPoint',[i_mnuSq,i_E0,i_B,i_N],...
            'Method', 'NonLinearLeastSquares',...
            'Algorithm','Levenberg-Marquardt','Display','off',...
            'Weights',(1./S.TBDIS));
        [fit1,gof,~] = fit(S.qU,S.TBDIS,myf,opts);
        par = [fit1.mb fit1.e0b fit1.bb fit1.nb 0 0];
        ci = confint(fit1,0.68);
        err = [(ci(2,1)-ci(1,1))/2 ...
            (ci(2,2)-ci(1,2))/2 ...
            (ci(2,3)-ci(1,3))/2 ...
            (ci(2,4)-ci(1,4))/2 ...
            0 0];
        chi2min = gof.sse;  % Sum of squares due to error - Equivalent to Chi2
end

% Definition of Amplitude / Background
mnuSq_Fit          = (S.mnuSq_i+par(1)); parA(1)=mnuSq_Fit;
mnuSq_FitError     = (err(1));

E0_Fit             = (par(2)); parA(2)=E0_Fit;
E0_FitError        = (err(2));

B_Fit              = (S.BKG_RateSec_i+par(3));  parA(3)=B_Fit;
B_FitError         = err(3);

N_Fit              = (par(4)); parA(4)=N_Fit;
N_FitError         = err(4);

parA(5)=0;
parA(6)=0;

switch display
    case 'ON'
        S.DisplayTDBInfo;
        fprintf(2,'-----------------------------------------------------------------------------------\n');
        fprintf('  Fit TBD Pixel \t= %g \n',PixelID);
        fprintf(2,'-----------------------------------------------------------------------------------\n');
        fprintf(2,'  m^2 \t = %g +- \t %g \t meV^2 \t - Init: %g \tmeV^2\n',mnuSq_Fit*1e6,mnuSq_FitError*1e6,i_mnuSq*1e6);
        fprintf(2,'  dE0 \t = %g +- \t %g \t meV \t - Init: %g \tmeV\n',E0_Fit*1e3,E0_FitError*1e3,i_E0*1e3);
        fprintf(2,'  B   \t = %g +- \t %g \t mcps \t - Init: %g \tmcps\n',B_Fit*1e3,B_FitError*1e3,i_B*1e3);
        fprintf(2,'  N   \t = %g +- \t %g \t \t - Init: %g \n', N_Fit, N_FitError,i_N);
        fprintf(2,'  Chi2\t = %g / %g dof \n',chi2min,S.nqU-4);
        fprintf(2,'-----------------------------------------------------------------------------------\n');
end

strtmp = sprintf('Fit and Residuals');
switch displayPlot
    case 'ON'
        fig = figure('Name',strtmp,'NumberTitle','off','rend','painters','pos',[10 10 800 600]);
    case 'OFF'
        fig = figure('Name',strtmp,'NumberTitle','off','rend','painters','pos',[10 10 800 600],'visible','off');
end
        subplot(2,2,[1 2])
        hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
             'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        %KurieCorr = Data(:,1).^2.*S.qUfrac;
%        hdata = errorbar(Data(:,1),sqrt(Data(:,2)./KurieCorr),sqrt(Data(:,3)./KurieCorr),...
%            'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
         hfit1 = plot(Data(:,1),tbd_modelint(par,S),...
             'Color','Black','LineWidth',1,'LineStyle','-');
%        hfit1 = plot(Data(:,1),sqrt(tbd_modelint(par,S)./KurieCorr),...
%            'Color','Black','LineWidth',1,'LineStyle','-');
 hfit3 = line([Data(1,1) Data(end,1)],[(S.BKG_RateSec_i+par(3))*S.TimeSec*S.qUfrac(end) (S.BKG_RateSec_i+par(3))*S.TimeSec*S.qUfrac(end)],'LineStyle','--','Color','Red');
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Counts/qU','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        mydata = sprintf('Data: m=%.2f eV - E0=%.2f eV \n',...
            sqrt(mnuSq_t),18575);
        myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n m=%.3f  +- %.3f eV \n E0=%.3f  +- %.3f meV \n B=%.3f  +- %.3f cps \n N=%.3f  +- %.3f',...
            chi2min,S.nqU-4,sqrt(abs(S.mnuSq_i+par(1))),err(1),(par(2))*1e3,err(2)*1e3,...
            (S.BKG_RateSec_i+par(3))*1e3,err(3)*1e3,par(4),err(4));
        mychi2 = sprintf('Data');

        l1 = legend([hdata  hfit1 hfit3],mychi2,myfit,'Background','Location','NorthEast') ;
        l1.FontSize = 11;
        axis([min(S.qU) max(S.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
        title(sprintf('Pixel %g - Data and Fit',S.FPD_Pixel),'FontSize',14);
        
        subplot(2,2,[3 4])
        hdata = errorbar(Data(:,1),Data(:,2)-tbd_modelint(par,S),Data(:,3),...
            'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
        hold on
        line([Data(1,1) Data(end,1)],[0 0],'LineStyle','--','Color','Blue');
        hold off
        grid on
        xlabel('qU (eV)','FontSize',14);
        ylabel('Residuals','FontSize',14);
        set(gca,'FontSize',14);
        set(gca,'yscale','lin');
        switch StatFluct
            case 'ON'
        axis([min(S.qU) max(S.qU)+1 min(Data(:,2)-tbd_modelint(par,S))*2 max(Data(:,2)-tbd_modelint(par,S))*2]);
        end
        
        if S.FPD_Pixel<10
            strname = sprintf('figures/%s_pixel00%g.eps',figprefix,S.FPD_Pixel);
        elseif S.FPD_Pixel<100
            strname = sprintf('figures/%s_pixel0%g.eps',figprefix,S.FPD_Pixel);
        else
            strname = sprintf('figures/%s_pixel%g.eps',figprefix,S.FPD_Pixel);
        end
        %publish_figure(1,strname);
        saveas(gcf,strname,'epsc');
        close

end
