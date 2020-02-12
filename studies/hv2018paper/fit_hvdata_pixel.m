function [par, err, chi2min, ndof] = fit_hvdata_pixel(varargin)
%
%
% Th. Lasserre - CEA Saclay
% Last Updated: January 2018
%

% Initialization
clear par ; clear parnomix;

% Parser
p = inputParser;
p.addParameter('mypub','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('KrObject','', @(x)isa(x,'Kr')); %Krypton Object
p.addParameter('CovMat', 'OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Pixel',1);
p.parse(varargin{:});

mypub             =    p.Results.mypub;
display           =    p.Results.display;
A                 =    p.Results.KrObject;
PixelID           =    p.Results.Pixel;
CovMat            =    p.Results.CovMat;

% Init
npar       = 5; i_E        = 0; i_W        = 0; i_Bkg      = 0; i_Phi0     = 0;i_N        = 0;
% Read Data File
[~, hvdatamap] = read_hvdata('Display','OFF','MK35',1972.4531);

% Read Data
Count    = (hvdatamap(PixelID,1,:));
CountErr = (hvdatamap(PixelID,2,:));
qUhv     = (hvdatamap(PixelID,3,:)); 
Data     = [qUhv(:),Count(:),CountErr(:)];
disp(Data);

% Init
A = init_hvdata_pixel('FPD_Pixel',PixelID);
A.L3_32_E_i = 30472.569 ; A.L3_32_W_i = 1.149;
LineE_i    = A.L3_32_E_i;
LineW_i    = A.L3_32_W_i;
LinePhi0_i = Count(1)-Count(end); A.L3_32_Phi0_i  = LinePhi0_i;
BKG_i      = Count(end);          A.BKG_RateSec_i = BKG_i;
A.qU       = squeeze(qUhv); A.SetKrDSBinning();
A.TimeSec  = 1; A.qUfrac  = ones(numel(A.qUfrac),1);
A.ComputeKrDS(); A.ComputeKrIS();

switch display
    case 'ON'
        %A.DisplayKrInfo;
        fprintf(2,'  Fit of Pixel %u with CovMat %s \n', A.FPD_Pixel, CovMat );
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  Initial Parameters \n');
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  Init E    \t= %.3f eV  \n',   LineE_i);
        fprintf(2,'  Init W    \t= %.3f eV \n',    LineW_i);
        fprintf(2,'  Init Phi0 \t= %.3f cps \n',   LinePhi0_i);
        fprintf(2,'  Init Bkg  \t= %.3f cps \n',   BKG_i );
        fprintf(2,'--------------------------------------------------------------\n');
end

switch CovMat
    case 'OFF'
        sigma = Data(:,3);
    case 'ON'
        StatErrMatrix  = (diag(Data(:,3))); %stat error only
        try 
            MyCovMat_tmp = load(sprintf('CovMatrix_hv2018paper_Pixel%u.mat', A.FPD_Pixel));
        catch
            fprintf('CovMat doesnt exists!! ' \n) 
        end
        MyCovMat = struct2cell(MyCovMat_tmp);
        sigma = MyCovMat{:,:} + StatErrMatrix.^2;
end

% Initializing Fit
% ndof = A.nqU-4;
% ParIni   = [i_E i_W i_Phi0 i_Bkg]; 
% Fitfunc = 'minimize.chi2';
% Theofunc = @A.FitCalculateKrIS;
% ParNames = {'E_Bias' 'W_Bias' 'Phi0_Bias' 'B_Bias'}; %Input for Theofunc
% ParNamesMinuit = cell2mat(strcat(ParNames,{' '}));   %For Minuit Instructions
% FitArg = {Data(:,2), sigma, Theofunc, ParNames };    %Cell: Counts, Uncertainties(Matrix), Theoriefunction, Parameterliste(strings) 
% MinuitArg = {'-n',ParNamesMinuit,'-c','set pri 1; min'};
% [par, err, chi2min, errmat] = fminuit(Fitfunc, ParIni, FitArg, MinuitArg{:});
% parValueCell=num2cell(par); 
% parArgIn={ParNames{:}; parValueCell{:}};
% fitresult = Theofunc(parArgIn{:}); %Counts

%Old Init without Fit class
ndof = A.nqU-4;
ParIni   = [i_E i_W i_Phi0 i_Bkg]; 
ParNames = {'E_Bias' 'W_Bias' 'Phi0_Bias' 'B_Bias'}; %Input for Theofunc
tmparg = sprintf(['set pri -1; min']);
DataL = {Data,A};
Args = {ParIni, DataL, '-n',ParNames,'-c ',tmparg};
[par, err, chi2min, errmat] = pminuit('chi2_hvdata_pixel',Args);
parValueCell=num2cell(par); 
parArgIn={ParNames{:}; parValueCell{:}};
Theofunc = @A.FitCalculateKrIS; fitresult = Theofunc(parArgIn{:}); %Counts

% Definition of Amplitude / Background
E_Fit             = (LineE_i+par(1)); parA(1)=E_Fit;
E_FitError        = (err(1));

W_Fit             = (LineW_i+par(2))*1e3; parA(2)=W_Fit;
W_FitError        = (err(2))*1e3;

Phi0_Fit          = (A.L3_32_Phi0_i+par(3))*(A.TimeSec.*A.qUfrac(1)); parA(3)=Phi0_Fit;
Phi0_FitError     = err(3)/(A.TimeSec.*A.qUfrac(1));

Offset_Fit        = (A.BKG_RateSec_i+par(4))*(A.TimeSec.*A.qUfrac(1));  parA(4)=Offset_Fit;
Offset_FitError   = err(4)/(A.TimeSec.*A.qUfrac(1));

parA(5) = 0;

switch display
    case 'ON'
        fprintf(2,'----------------------------------------------------\n');
        fprintf(2,'Fit Kr83m L3-32 Line \n');
        fprintf(2,'----------------------------------------------------\n');
        fprintf(2,'  E \t= %.3f \t +- \t %g \t eV \n',E_Fit,E_FitError);
        fprintf(2,'  W \t= %.3f \t +- \t %g \t meV \n',W_Fit,W_FitError);
        fprintf(2,'  Phi0 \t= %.3f +- \t %g \t \n',Phi0_Fit,Phi0_FitError);
        fprintf(2,'  Bkg \t= %.3f \t +- \t %g \t cps \n',Offset_Fit,Offset_FitError);
        %fprintf(2,'  N \t= %g \t +- \t %g \t fixed \n',(par(5)),err(5));
        fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
        fprintf(2,'----------------------------------------------------\n');
end

%% Plot Results
str = sprintf('Fit and Residuals - PixelID %.0f',PixelID);
fig = figure('Name',str,'NumberTitle','off','rend','painters','pos',[10 10 1400 800]);

subplot(2,2,1)
hdata = errorbar(Data(:,1),Data(:,2),Data(:,3),...
    'ks','MarkerSize',3,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on

%hfit1 = plot(Data(:,1),model_hvdata_pixel(par,A),...
hfit1 = plot(Data(:,1),fitresult,...
    'Color','Black','LineWidth',1,'LineStyle','-');
hfit3 = line([Data(1,1) Data(end,1)],[Offset_Fit Offset_Fit],'LineStyle','--','Color','Red');
hold off
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Counts/qU','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
mydata = sprintf('Data Pixel %.0f: E=%.2f eV - W=%.2f eV \n',...
    PixelID,LineE_i,LineW_i);
myfit = sprintf(' Fit: \\chi2 / dof=%.1f/%.0f \n E=%.3f \t\\pm %.3f eV \n W=%.0f \t\\pm %.0f meV \n A=%.3f \t\\pm %.3f \n O=%.3f \t\\pm %.3f',...
    chi2min,A.nqU-4,E_Fit,E_FitError,W_Fit,W_FitError,...
    Phi0_Fit,Phi0_FitError,Offset_Fit,Offset_FitError);
mychi2 = sprintf('Data');
legend([hdata  hfit1 hfit3],mychi2,myfit,'Offset','Location','NorthEast') ;
axis([min(A.qU) max(A.qU)+1 0.*min(Data(:,2)) max(Data(:,2))*1.2])
title(sprintf('KATRIN Kr83m-g - {%s} Pixel %.0f',...
    A.TD,PixelID),'FontSize',14);

subplot(2,2,3)
hdata = errorbar(Data(:,1),Data(:,2)-fitresult,Data(:,3),...
    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Blue');
hold off
grid on
xlabel('qU (eV)','FontSize',14);
ylabel('Residuals','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
axis([min(A.qU) max(A.qU)+1 min(Data(:,2)-fitresult)*2 max(Data(:,2)-fitresult)*2]);

subplot(2,2,2)
hfit2 = plot(A.Te,A.KrDS,...
    'Color','Black','LineWidth',1,'LineStyle','-');
grid on
xlabel('Te (eV)','FontSize',14);
ylabel('Counts/Te','FontSize',14);
set(gca,'FontSize',14);
set(gca,'yscale','lin');
title(sprintf('Differential Spectrum - Doppler %s (fixed \\sigma=%.1f eV)',...
    A.DopplerEffectFlag,A.DE_sigma),'FontSize',14);

subplot(2,2,4)
[~, cor] = corplot(errmat(1:4,1:4));
xticklabels({'Position','Width','Amplitude','Offset'})
yticklabels({'Position','Width','Amplitude','Offset'})
disp(cor);

str = sprintf('Run33149-PixelID%.0f.png',PixelID);
saveas(gcf,str)

end
