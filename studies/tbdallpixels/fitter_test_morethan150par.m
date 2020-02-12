clear par ; clear parnomix;
tic; varargin = {};
fitter = 'mat';
% Parser
p = inputParser;
p.addParameter('Run',110,@(x)isfloat(x));
p.addParameter('MyMode','Simulation',@(x)ismember(x,{'Bootcamp','Simulation'}));
p.addParameter('pub','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('display','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('StatFluct','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('nPixels',3,@(x)isfloat(x) && x<150);
p.addParameter('FixNorms','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('Backgrounds','READ',@(x)ismember(x,{'FIT','READ'}));
p.addParameter('mnuSq_t',0.,@(x)isfloat(x));
p.addParameter('FPDView','OFF',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});
Run         =    p.Results.Run;
MyMode      =    p.Results.MyMode;
pub         =    p.Results.pub;
display     =    p.Results.display;
StatFluct   =    p.Results.StatFluct;
nPixels     =    p.Results.nPixels;
FixNorms    =    p.Results.FixNorms;
Backgrounds =    p.Results.Backgrounds;
mnuSq_t     =    p.Results.mnuSq_t;
FPDView     =    p.Results.FPDView;

switch MyMode
    case 'Bootcamp'
        switch Run
            case 110
                A = init_bootcampR110('PixelID',0,'mnuSq_i',mnuSq_t,'nPixels',nPixels);
                D = load('data/runmat1_10.mat');
                figprefix = 'bc_fitallsim_run1_10';
        end
end

switch MyMode
    case 'Bootcamp'
        quIndex        = 1:3:A.nPixels*3;
        countIndex     = quIndex + 1;
        coutErrIndex   = quIndex + 2;
        A.qUpixel(1:A.nPixels,:) = D.TBDISallPixels(:,quIndex(1:A.nPixels))'; 
        switch display
            case 'ON'
                fprintf(2,'----------------------------------------\n');
                fprintf(2,'qU Matrix - %g Pixels\n',A.nPixels);
                disp(num2str(A.qUpixel(1:A.nPixels,:)));
                fprintf(2,'----------------------------------------\n');
        end
    case 'Simulation'
        A.qUpixel(1:A.nPixels,:) = repmat(A.qU',A.nPixels,1);
end

% Init Normalization
fitr = load(sprintf('results/bc_fitallpixels_run%g_G.mat',Run));

% Init
i_mnuSq    = A.mnuSq_i;%(0.2*randn)^2;
i_E0       = 0;%0.05*randn;
i_N        = fitr.ParnPixels(1:A.nPixels,5)';
i_B        = fitr.ParnPixels(1:A.nPixels,4)';

switch MyMode
    case 'Bootcamp'
        Data = [...
            reshape(D.TBDISallPixels(:,quIndex(1:A.nPixels)),1,A.nqU*A.nPixels)',...
            reshape(D.TBDISallPixels(:,countIndex(1:A.nPixels)),1,A.nqU*A.nPixels)',...
            reshape(D.TBDISallPixels(:,coutErrIndex(1:A.nPixels)),1,A.nqU*A.nPixels)'];
        DataL = {Data,A};
            case 'Simulation'
        A.ComputeTBDDSallPixels(); A.ComputeTBDISallPixels();
        switch StatFluct
            case 'ON'
                AddStatFluctTBDISallPixels(A);
        end
        Data = [...
            reshape(repmat(A.qU,1,A.nPixels),1,A.nqU*A.nPixels)' , ...
            reshape((A.TBDISallPixels),1,A.nqU*A.nPixels)', ...
            reshape(sqrt((A.TBDISallPixels)),1,A.nqU*A.nPixels)'];
        DataL = {Data,A};
end

ParIni    = [i_mnuSq i_E0 i_N i_B];

% Build Name Parameter List
parnames = ['mSq E0'];
parnamescor = [sprintf('''mSq'',''E0''')];
for k=1:1:A.nPixels
tmp         = sprintf(' N%g',k);
parnames    = [parnames tmp];
tmpcor      = sprintf(',''N%g''',k);
parnamescor = [parnamescor tmpcor];
end

for k=1:1:A.nPixels
    tmp         = sprintf(' B%g',k);
    parnames    = [parnames tmp];
    tmpcor      = sprintf(',''B%g''',k);
    parnamescor = [parnamescor tmpcor];
end

%% MINIMIZATION HERE XXX
options = optimoptions('fminunc','Algorithm','quasi-newton',...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-7,...
    'MaxFunctionEvaluations',1e6,'FiniteDifferenceType','central');
TBDfun = @(xx) chi2_tbd_allpixelssim_m(xx,DataL);
[par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,ParIni,options);
if exitflag ~= 1
    [par,chi2min,exitflag,output,grad,Hessian] = fminunc(TBDfun,par,options);
end
errmat = 0.5*Hessian;
varcov = inv(errmat);
err = sqrt(diag(varcov));


%% Save Results
strtmp = sprintf('./results/bc_fitallpixels_run%g_sim_%gp.mat',Run,nPixels);
save(strtmp,'par','err','chi2min','errmat');

%% Plot Results
figure(1)
subplot(2,1,1)
hdata = errorbar(Data(1:A.nqU*A.nPixels,1),Data(1:A.nqU*A.nPixels,2),Data(1:A.nqU*A.nPixels,3),...
    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
model = model_tbdallpixelssim_m(par,A);

for (i=0:1:A.nPixels-1)
    hfit1 = plot(Data(i*A.nqU+1:i*A.nqU+A.nqU,1),model(i*A.nqU+1:i*A.nqU+A.nqU),...
        'Color','Red','LineWidth',1,'LineStyle','-');
end

hold off
grid on
xlabel('qU (eV)','FontSize',10);
ylabel('Counts','FontSize',10);
title(sprintf('KATRIN TBD - %s - Pixel 1 to %.0f - Doppler %s',A.TD,A.nPixels,A.DopplerEffectFlag));
set(gca,'FontSize',12);
set(gca,'yscale','lin');
mydata = sprintf('Data: NuMassSq=%.4f eV^2 - E0=%.4f meV \n',A.mnuSq_i,A.Q);
myfit = sprintf('Fit: mSq=%.4f \\pm %.4f eV - E0=%.4f \\pm %.4f eV',(i_mnuSq+par(1)),err(1),i_E0+par(2),err(2));
mychi2 = sprintf('\\chi2 / dof=%.1f/%.0f\n',chi2min,A.nqU*A.nPixels-2);
legend([hdata  hfit1],mychi2,myfit,'Location','NorthEast') ;% legend(a,'boxoff');
axis([min(A.qU) max(A.qU)+1 0.7*min(Data(:,2)) max(Data(:,2))*1.2])

subplot(2,1,2)
hdata = errorbar(Data(1:A.nqU*A.nPixels,1),Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels),Data(1:A.nqU*A.nPixels,3),...
    'ks','MarkerSize',5,'MarkerFaceColor',.9*[1 1 1],'LineWidth',1);
hold on
line([Data(1,1) Data(end,1)],[0 0 ],'LineStyle','--','Color','Red');
hold off
grid on
xlabel('qU (eV)','FontSize',10);
ylabel('Residuals','FontSize',10);
set(gca,'FontSize',12);
set(gca,'yscale','lin');
%axis([min(A.qU) max(A.qU)+1 min(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))*2 max(Data(1:A.nqU*A.nPixels,2)-model(1:A.nqU*A.nPixels))*2])
%publish_figure(1,'figures/fit_tbdallpixelssim1.eps');

figure(2)
title(sprintf('Correlation Matrix of %g Fitted Parameters',npar));
        [h cor] = corplot(errmat(1:numel(ParIni),1:numel(ParIni)));
%         xticklabels({'mSq','E0','N1','N2','N3'})
%         yticklabels({'mSq','E0','N1','N2','N3'})
%        disp(parnamescor);
%publish_figure(2,'figures/fit_tbdallpixelssim2.eps');

figure(3)
title(sprintf('TBD Integral Spectra and Fit -  %g Pixels',A.nPixels));
hold on
[x, y] = meshgrid(A.qU,1:1:A.nPixels);
Rfit  = reshape(model(1:A.nqU*A.nPixels),A.nqU,A.nPixels);
Rdata = reshape(Data(1:A.nqU*A.nPixels,2),A.nqU,A.nPixels);
ribbon(Rfit);colormap('summer')
hold on
stem3(Rdata,'Marker','s',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'LineStyle','none');
hold off
ylabel('qU')
xlabel('Pixel')
zlabel('Integral Spectrum Fit')
hold off
view([-100 45])
%publish_figure(3,'figures/fit_tbdallpixelssim3.eps');

figure(4)
title(sprintf('TBD Integral Spectra and Fit Residuals -  %g Pixels',A.nPixels));
hold on
[x, y] = meshgrid(A.qU,1:1:A.nPixels);
Rfit     = reshape(model(1:A.nqU*A.nPixels),A.nqU,A.nPixels);
Rdata    = reshape(Data(1:A.nqU*A.nPixels,2),A.nqU,A.nPixels);
RdataErr = reshape(Data(1:A.nqU*A.nPixels,3),A.nqU,A.nPixels);
ribbon(Rfit-Rfit);colormap('summer')
hold on
stem3((Rdata-Rfit)./RdataErr,'Marker','s',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'LineStyle','none');
hold off
ylabel('qU')
xlabel('Pixel')
zlabel('Integral Spectrum Fit Residuals')
hold off
view([180 0])
%publish_figure(4,'figures/fit_tbdallpixelssim4.eps');
