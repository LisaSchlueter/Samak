% test implementation of penning trap background slope for knm2
% lisa, jan '21
Plot = 'ON';
DataType = 'MC';
range = 40;
BKG_PtSlope =1.7*1e-06; %3e-06;%1.7*1e-06; % background rate increase over time, in cps per second 
freePar = 'mNu E0 Bkg Norm BkgPTSlope';

savedir = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];
savename = sprintf('%sknm2_PngBkg_%s_%.0feV_Bpng-%.1fmucpsPers.mat',...
    savedir,DataType,range,BKG_PtSlope*1e6);

if exist(savename,'file')
    load(savename);
else
    %% settings
    SigmaSq =  0.0124+0.0025;
    RunAnaArg = {'RunList','KNM2_Prompt',...
        'DataType','Real',...
        'fixPar',freePar,...
        'minuitOpt','min ; minos',...
        'FSDFlag','KNM2',...
        'ELossFlag','KatrinT2A20',...
        'SysBudget',38,...
        'AnaFlag','StackPixel',...
        'chi2','chi2Stat',...
        'TwinBias_Q',18573.7,...
        'NonPoissonScaleFactor',1,...
        'FSD_Sigma',sqrt(SigmaSq),...
        'TwinBias_FSDSigma',sqrt(SigmaSq),...
        'RingMerge','Full',...
        'PullFlag',99,...
        'RadiativeFlag','ON'};
    
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    if ismember(A.DataType,{'Twin','MC'})
        A.ModelObj.RFBinStep = 0.01;
        A.ModelObj.InitializeRF;
    end
    
    A.ModelObj.nRuns = A.nRuns;
    %% model w/o time dependent bkg slope from penning track
    % reference fit w/o bkg slope from penning trap in model
    A.ModelObj.BKG_PtSlope = 0; % in counts per second (1.7 ± 3.0) µcps/s
    A.InitModelObj_Norm_BKG;
    
    A.ModelObj.mnuSq_i = 0.1;
    A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
    TBDIS_01eV2 = A.ModelObj.TBDIS;
    
    A.ModelObj.mnuSq_i = -0.1;
    A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
    TBDIS_m01eV2 = A.ModelObj.TBDIS;
    
    A.ModelObj.mnuSq_i = 0 ;
    A.ModelObj.ComputeTBDDS; A.ModelObj.ComputeTBDIS;
    TBDIS_i = A.ModelObj.TBDIS;

    if strcmp(DataType,'MC')
        A.RunData.TBDIS = TBDIS_i;
    end
    
    A.Fit;
    FitResult_i = A.FitResult;
    
    %% model with time dependent bkg slope from penning track
    A.ModelObj.BKG_PtSlope = BKG_PtSlope;
    
    A.ModelObj.ComputeTBDIS;
    TBDIS_tot = A.ModelObj.TBDIS;
    TBDIS_BkgPngSlope = TBDIS_tot - TBDIS_i;
    
    A.Fit;
    FitResult_png = A.FitResult;
    
    %% prep for display
    
    qU = A.ModelObj.qU-18574;
    Time = A.ModelObj.TimeSec.*A.ModelObj.qUfrac;
    BkgRate = A.ModelObj.BKG_RateSec;
    nqU = A.ModelObj.nqU;
    exclDataStart = A.exclDataStart;
    Q_i = A.ModelObj.Q_i;
    MakeDir(savedir);
    save(savename,'qU','Time','BkgRate','nqU','TBDIS_tot','TBDIS_i',...
        'TBDIS_BkgPngSlope','exclDataStart','RunAnaArg','FitResult_i','FitResult_png','Q_i',...
        'TBDIS_m01eV2','TBDIS_01eV2');
end

Mode = 'Rate';
switch Mode
    case 'Counts'
        ytot =  TBDIS_tot;
        y_i  =  BkgRate.*Time;% TBDIS_i;
        y_png = TBDIS_BkgPngSlope;
        yStr = 'Counts';
    case 'Rate'
        
        ytot = TBDIS_tot./Time;
        y_i  =  repmat(BkgRate,nqU,1);
        y_png = TBDIS_BkgPngSlope./Time;
        yStr = 'Rate (cps)';
end

%% Plot 1 : spectrum
if strcmp(Plot,'ON')
    f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.8]);
    s1 = subplot(4,1,[1,2]);
    ptot =plot(qU,ytot,'-.','LineWidth',2.5,'Color',rgb('Gray'));
    hold on;
    pi= plot(qU,y_i,':','LineWidth',2,'Color',rgb('LightGreen'));
    ppng =plot(qU,y_png,'-','LineWidth',2,'Color',rgb('DodgerBlue'));
    PrettyFigureFormat;
    set(gca,'YScale','log');
    ylabel(yStr);
    ylim([min(y_png(exclDataStart:end))-1e-04,max(ytot(exclDataStart:end))+1])
    leg = legend([ptot,pi,ppng],...
        sprintf('\\beta + {\\it B}_{const.} + {\\it B}_{pt.}'),...
        sprintf('{\\it B}_{const.}'),...
        sprintf('{\\itB}_{pt.}: \\alpha = %.1f \\mucps/s',BKG_PtSlope*1e6));
    leg.EdgeColor = rgb('LightGray');
    xticklabels('')
    ax1 = gca;
    t = title(sprintf('Penning trap induced background rate increase {\\itB}_{pt.} = %.1e cps/s \n%s:  \\Delta{\\itm}_\\nu^2 = %.3f eV^2 , \\Delta{\\itE}_0 = %.3f eV',...
        BKG_PtSlope,DataType,FitResult_png.par(1)-FitResult_i.par(1),FitResult_png.par(2)-FitResult_i.par(2)),'FontWeight','normal','FontSize',get(ax1,'FontSize'));
    ax1.Position = [ax1.Position(1) ax1.Position(2)-0.04, ax1.Position(3),ax1.Position(4)];
 
    
    s2 = subplot(4,1,3);
    ppng = plot(qU,TBDIS_tot./TBDIS_i,'-','LineWidth',2,'Color',rgb('DodgerBlue'));%.*TBDIS_i(exclDataStart)./TBDIS_tot(exclDataStart),'-');
    hold on;
    p05eV  = plot(qU,TBDIS_01eV2./TBDIS_i,'-.','LineWidth',2,'Color',rgb('Orange'));%.*TBDIS_i(exclDataStart)./TBDIS_1eV2(exclDataStart),':');
    pm05eV = plot(qU,TBDIS_m01eV2./TBDIS_i,':','LineWidth',2,'Color',rgb('FireBrick'));
    PrettyFigureFormat;
    ylabel(sprintf('Ratio'))
    leg = legend([ppng,p05eV,pm05eV],sprintf('{\\itB}_{pt.}: \\alpha = %.1f \\mucps/s',BKG_PtSlope*1e6),...
        sprintf('{\\itm}_\\nu^2 = 0.1 eV^2'),...
        sprintf('{\\itm}_\\nu^2 = -0.1 eV^2'));
    ax1 = gca;
    ax1.Position = [ax1.Position(1) ax1.Position(2)-0.04, ax1.Position(3),0.04+ax1.Position(4)];
    xticklabels('')
    yticks([0.999,1,1.001]);
    leg.EdgeColor = rgb('LightGray');
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.9]));
    
    s3 = subplot(4,1,4);
    b1 = bar(qU,Time./(60*60),'FaceColor',rgb('DodgerBlue'),'EdgeColor',rgb('DodgerBlue'));
    PrettyFigureFormat;
    xlabel('Retarding potential -18574 (eV)');
    linkaxes([s1,s2,s3],'x');
    xlim([-range,137])
    ylabel('Time (h)');
    ax1 = gca;
    ax1.Position = [ax1.Position(1) ax1.Position(2)-0.04, ax1.Position(3),0.04+ax1.Position(4)];
    
    plotdir = strrep(savedir,'results','plots');
    MakeDir(plotdir);
    plotname =  sprintf('%sknm2_PngBkg_%s_%.0feV_Bpng-%.1fmucpsPers.png',...
        plotdir,DataType,range,BKG_PtSlope*1e6);
    print(gcf,plotname,'-dpng','-r300');
end
%% nu-mass bias
fprintf('reference fit: mnu^2 = %.3f eV^2 +-  %.3f eV^2, E0 = %.3f eV \n',...
    FitResult_i.par(1),FitResult_i.err(1),FitResult_i.par(2));
fprintf('BslopePng fit: mnu^2 = %.3f eV^2 +- %.3f eV^2 , E0 = %.3f eV \n',...
    FitResult_png.par(1),FitResult_i.err(1),FitResult_png.par(2));


