% test conv function by convolving known functions -> sanity check!
Mode = 'Step';%'Gauss'; % 'Step'
BinStep = 0.01;
savedir = [getenv('SamakPath'),'knm2ana/knm2_ResponseFunction/results/'];
savename = sprintf('%sTestConv_%s_%.2fStep.mat',savedir,Mode,BinStep);

if exist(savename,'file') 
    load(savename)
else
    e = -50:BinStep:50;
    Gaussian = @(x,mu,sigma) 1./sqrt(2*pi*sigma^2).*exp(-0.5*(x-mu).^2./sigma.^2);
    StepFun = @(e) 0+(e>=6.5 & e<=7)./2;
    
    fun1 = Gaussian(e,5,1);
    switch Mode
        case 'Step'
            fun2 = StepFun(e);
            legStr = sprintf('Step func.:   \\mu = 6.5    , width   = 1');
            f2 = @(e) StepFun(e);
        case 'Gauss'
            fun2 = Gaussian(e,8,2);
            legStr = sprintf('Gaussian 2:   \\mu = %.1f   , \\sigma^2  = %.1f',8,2^2);
            f2 = @(e) Gaussian(e,8,2);
    end
    tic;
    Convfun = conv(fun1,fun2,'same').*(e(2)-e(1));
    toc;
    
    FitRes =fit(e',Convfun','gauss1');
    SigmaConv = FitRes.c1/sqrt(2);
    
    %% new convolution 1) --> accurate, but to slow for response function :(
    f1 = @(e) Gaussian(e,5,1)';
    k = e';
    Intfun = @(x) f1(x).*f2(k-x);
    tic;
    ConvfunInt = integral(Intfun,0,100,'ArrayValued',1);
    toc;
    %% new convolution 2)
    Eint = linspace(0,max(e),1e4); % for integration
    tic;
    ConvfunInt2 = IntConv(e,Eint,f1,f2);
    toc;
    %% save
    MakeDir(savedir);
    save(savename,'e','fun1','fun2','Convfun','ConvfunInt','ConvfunInt2','SigmaConv','FitRes','legStr','f1','f2');
end
%% sanity plot
figure('Units','normalized','Position',[0.1,0.1,0.45,0.4])
plot(e,fun1,e,fun2,e,Convfun,'LineWidth',2);
hold on;
plot(e,ConvfunInt,'-.','LineWidth',2);
plot(e,ConvfunInt2,':','LineWidth',2);
leg = legend(sprintf('Gaussian:     \\mu = %.1f   , \\sigma^2  = %.1f',5,1),...
    legStr,...
    sprintf('Convolution {\\itconv}'),...%  \\mu = %.1f , \\sigma^2  = %.1f',FitRes.b1,SigmaConv^2),...
    sprintf('Convolution {\\itintegral}'),...% max. diff. = %.2g',max(abs(ConvfunInt-Convfun'))),...
    sprintf('Convolution {\\itsimpsons}'),...%: max. diff. = %.2g',max(abs(ConvfunInt2-Convfun'))),...
    'EdgeColor',rgb('Silver'));
PrettyFigureFormat;
xlim([0 20]);
ylim([0 max([fun2,fun1])*1.1]);
title(sprintf('Bin width for {\\itconv}() = %.2f',BinStep),'FontWeight','normal');

plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(gcf,plotname);
fprintf('save plot to %s \n',plotname);