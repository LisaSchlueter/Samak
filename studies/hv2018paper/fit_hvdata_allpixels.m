function [par,err,chi2min,ndof] = fit_hvdata_allpixels(varargin)
%
%            FIT INTEGRAL SPECTRUM
%              Fit a 83mKr Lines
%           All Pixels Simultaneously 
%             For up to 73 Pixels
%  
%         Common E0 and W for all Pixels
%          Individual Phi0 and Offset
%
%          Th. Lasserre - CEA Saclay
%                October 2017
%                       &
%            L.Schlueter - MPP Muenchen
%                Feburary 2018

% Initialization
clear par ; clear parnomix;

% Parser
p = inputParser;
p.addParameter('fign',1,@(x)isfloat(x) && x>0);
p.addParameter('pub',@(x)ismember(x,{'ON','OFF','YES'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('CovMat', 'OFF', @(x)ismember(x,{'ON','OFF'}));
p.addParameter('fixPar', 'none', @(x)ismember(x,{'none','Phi0', 'BKG'})); 
p.parse(varargin{:});
fign       =    p.Results.fign;
pub        =    p.Results.pub;
display    =    p.Results.display;
CovMat     =    p.Results.CovMat;
fixPar     =    p.Results.fixPar;

% Parametrization: True Value
A = init_hvdata_pixel('FPD_Pixel',0, 'nPixels',40);

A.TimeSec = 1;
A.qUfrac  = ones(numel(A.qUfrac),1);
%A.SetKrDSBinning(A.qU(1),A.qU(end),1e-2);
A.DisplayKrInfo;

% Read Data File - old style, 1 qU vector for all pixels
[~, hvdatamap] = read_hvdata('Display','OFF','MK35',1972.4531);
% Read Data File - new style, 1 qU vector per pixel
Data = [reshape(permute(hvdatamap(1:A.nPixels,3,1:end),[3 2 1]),[],1),...
    reshape(permute(hvdatamap(1:A.nPixels,1,1:end),[3 2 1]),[],1),...
    1.*reshape(permute(hvdatamap(1:A.nPixels,2,1:end),[3 2 1]),[],1)];

StatErrMatrix  = (diag(Data(:,3))); %stat error only

switch CovMat
    case 'OFF'
        FitCovErrMatrix = StatErrMatrix.^2;
    case 'ON'
        try
        MultiCM_tmp = load(sprintf('MultiCovMatrix_hv2018paper_all%uPixels.mat',A.nPixels));
        MultiCM = struct2cell(MultiCM_tmp);
        catch
            fprintf('Multi CovMat doesnt exists!! ' \n)
        end
         FitCovErrMatrix = MultiCM{:,:} + StatErrMatrix.^2; 
end
%DataL = {Data,A, FitCovErrMatrix};

%Init
A.L3_32_Phi0allPixels_i  = (hvdatamap(1:A.nPixels,1,1)-hvdatamap(1:A.nPixels,1,end))';
A.BKG_RateSecallPixels_i = (hvdatamap(1:A.nPixels,1,end))';
A.SetLineInitParallPixels;
A.L3_32_E_i = 30473; 
A.L3_32_W_i = 1.15;

A.ComputeKrDSallPixels();
A.ComputeKrISallPixels();

LineE_i    = A.L3_32_E_i;
LineW_i    = A.L3_32_W_i;
LineBKG_i  = A.BKG_RateSecallPixels_i;
LinePhi0_i = A.L3_32_Phi0allPixels_i;

switch display
    case 'ON'
        fprintf(2,'----- FIT INIT ------\n');
        fprintf(2,'CovMat %s \n',CovMat);
        fprintf(2,'E_i = %.3f eV \n', LineE_i);
        fprintf(2,'W_i = %.3f eV \n', LineW_i);
        for i=1:A.nPixels
            fprintf(2,'Pixel %u \n', i);
            fprintf(2,'Phi_i= %.3f  counts \n', LinePhi0_i(i));
            fprintf(2,'BKG_i = %.3f counts \n', LineBKG_i(i));
        end
end

[ParIni fix npar] = A.FitInit(fixPar); %Init
ParNamesMinuit = A.BuildFitStrings; %Build string with all Parameter
ParNames = {'E_Bias', 'W_Bias', 'Phi0_Bias', 'B_Bias'}; %for Model input

% Fit New
Fitfunc = 'minimize.chi2';
Theofunc = @A.FitCalculateKrISallPixels;
FitArg = {Data(:,2), FitCovErrMatrix, Theofunc, ParNames }; %Cell: Counts, Uncertainties(Matrix), Theoriefunction, Parameterliste(strings) 
tmpArgMinuit = sprintf('set pri 2; %s ; min', fix);
MinuitArg = {'-n',ParNamesMinuit,'-c',tmpArgMinuit};
fprintf(2,tmpArgMinuit);
[par, err, chi2min, errmat] = fminuit(Fitfunc, ParIni, FitArg, MinuitArg{:});

% Fit Old
% tmparg = sprintf(['set pri 1; min']);
% DataL = {Data,A};
% Args = {ParIni, DataL, '-n',ParNames,'-c ',tmparg};
% [par, err, chi2min, errmat] = pminuit('chi2_hvdata_allpixels',Args);

% Display
switch display
    case 'ON'
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf('  Processing ....\n');
        fprintf(2,'--------------------------------------------------------------\n');
        fprintf(2,'  E \t= %.3f \t +- \t %g \t eV \n',  (LineE_i+par(1)),err(1));
        fprintf(2,'  W \t= %.3f \t +- \t %g \t eV \n',  (LineW_i+par(2)),err(2));
        for i=1:1:A.nPixels
        fprintf(2,' Pixel %u Phi0 \t= %.3f \t +- \t %g \t eV \n', i, (LinePhi0_i(i)+par(2+i)),err(2+i));
        fprintf(2,' Pixel %u Off. \t= %.3f \t +- \t %g \t eV \n', i, (LineBKG_i(i)+par(A.nPixels+2+i)),err(A.nPixels+2+i));
        end  
        %fprintf(2,'  N \t= %.3f \t +- \t %g \t  \n',    par(5),err(5));
        ndof = A.nqU*A.nPixels-(A.nPixels*2+2);
        fprintf(2,'  Chi2 \t= %g / %g dof \n',chi2min,ndof);
        fprintf(2,'--------------------------------------------------------------\n');
end

%save results
switch CovMat
    case 'OFF'
        str1 = sprintf('./results/FitResults_hvdata_all_%s_newB-FieldCorr',date);
    case 'ON'
        str1 = sprintf('./results/FitResultsCM_hvdata_all_%s',date);
end
str2 = num2str(A.nPixels); str3 = 'Pixels.mat';
filename = strcat(str1, str2, str3);
save(filename);
%% Plot Results
orange=1/255*[255,140,0];
figure(fign+999)
subplot(2,1,1)
hdata = errorbar(Data(1:A.nqU*A.nPixels,1),Data(1:A.nqU*A.nPixels,2),Data(1:A.nqU*A.nPixels,3),...
    's', 'MarkerSize',5,'LineWidth',1);
hold on
%     hfit1 = plot(A.Te,KrLineModelDiff4par(par)./trapz(A.Te,KrLineModelDiff4par(par)).*trapz(A.qU,KrLineModel4par(par))/4,...
%         'LineWidth',1,'LineStyle','--','Color','Black');
model = model_hvdata_allpixels(par,A);

for (i=0:1:A.nPixels-1)
    hfit1 = plot(Data(i*A.nqU+1:i*A.nqU+A.nqU,1),model(i*A.nqU+1:i*A.nqU+A.nqU),...
        'Color',orange,'LineWidth',1,'LineStyle','-');
end

hold off
grid on
xlabel('qU (eV)','FontSize',18);
ylabel('Counts','FontSize',18);
title(sprintf('KATRIN Krypton gas - %s - 40 Pixel Fit (CM)', A.TD), 'FontSize', 16);
%set(gca,'FontSize',12);
set(gca,'yscale','lin');
mydata = sprintf('Data: E=%.2f eV - W=%.2f eV \n',LineE_i,LineW_i);
myfit = sprintf('E=%.3f \\pm %.3f eV \nW=%.3f \\pm %.4f eV',LineE_i+par(1),err(1),LineW_i+par(2),err(2));
mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,ndof);
leg1 = legend([hdata  hfit1],mychi2,myfit,'Location','NorthEast') ;% legend(a,'boxoff');
leg1.FontSize = 14;
axis([min(A.qU) max(A.qU) 0.7*min(Data(:,2)) max(Data(:,2))*1.2])

subplot(2,1,2)
%hdata = errorbar(Data(1:A.nqU*A.nPixels,1),(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels)),Data(1:A.nqU*A.nPixels,3),...
%    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hdata = plot(Data(1:A.nqU*A.nPixels,1),(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))./Data(1:A.nqU*A.nPixels,3),...
      's','MarkerSize',3,'LineWidth',1);
hold on
line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','-','Color','Red', 'LineWidth',2);

hold off
grid on
xlabel('qU (eV)','FontSize',18);
ylabel('norm. Residuals','FontSize',17);
%set(gca,'FontSize',12);
set(gca,'yscale','lin');
%axis([min(A.qU) max(A.qU)+1 min(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))*2 max(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))*2])
xlim([min(A.qU) max(A.qU)]);
export_fig './figures/Fit40PixelsCM_spectrum.BCorrfig'
export_fig './figures/Fit40PixelsCM_spectrum_BCorr.pdf'


% % Create a grid of x and y points
% [x, y] = meshgrid(A.qU,1:1:A.nPixels);
% Rfit  = reshape(model(1:A.nqU*A.nPixels),A.nqU,A.nPixels);
% Rdata = reshape(Data(1:A.nqU*A.nPixels,2),A.nqU,A.nPixels);
% figure;
% ribbon(Rfit);colormap('summer')
% hold on
% stem3(Rdata,'Marker','s',...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','k',...
%     'LineStyle','none');
% hold off
% ylabel('qU')
% xlabel('Pixel')
% zlabel('Integral Spectrum Fit')
end
