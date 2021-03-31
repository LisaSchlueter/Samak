function relic_global(varargin)

    %% Settings
    p=inputParser;
    p.addParameter('eta',2.8e9,@(x)isfloat(x));
    p.addParameter('Params','TDR',@(x)ismember(x,{'TDR','Formaggio','KNM1','KNM2'}));
    p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
    p.addParameter('Init_Opt','',@(x)iscell(x) || isempty(x));
    p.addParameter('E0',18575,@(x)isfloat(x));                                         % Endpoint in eV
    p.addParameter('Syst','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('range',40,@(x)isfloat(x));
    p.addParameter('saveplot','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('annotations','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    eta         = p.Results.eta;
    Params      = p.Results.Params;
    fitPar      = p.Results.fitPar;
    Init_Opt    = p.Results.Init_Opt;
    E0          = p.Results.E0;
    Syst        = p.Results.Syst;
    range       = p.Results.range;
    saveplot    = p.Results.saveplot;
    annotations = p.Results.annotations;
    
    if strcmp(Syst,'ON')
        Chi2opt='chi2CMShape';
    else
        Chi2opt='chi2Stat';
    end
    
    switch Params
        case 'TDR'
            initfile=@ref_RelicNuBkg_DesignReport;
        case 'Formaggio'
            initfile=@ref_RelicNuBkg_Formaggio;
        case 'KNM1'
            initfile=@ref_RelicNuBkg_KNM1;
    end

    %% Data
    D = RunAnalysis('RunNr',100,...
        'RecomputeFakeRun','ON',...
        'Init_Opt',[Init_Opt,{'eta_i',eta}],...
        'FakeInitFile',initfile,...
        'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'TwinBias_Q',E0,...
        'RingList',1:14,...
        'fixPar',fitPar,...                   % free Parameter!!
        'NonPoissonScaleFactor',1,...         % background uncertainty are enhanced
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','Sibille0p5eV',...          % final state distribution                        !!check ob initfile hier überschrieben wird
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','FSD',...
        'SynchrotronFlag','OFF',...
        'AngularTFFlag','OFF');
    
    %D.InitModelObj_Norm_BKG('RecomputeFlag','ON');

    D.exclDataStart = D.GetexclDataStart(range);
    %% No relics

    if isempty(Init_Opt)
        R = RunAnalysis('RunNr',1,...
            'FakeInitFile',initfile,...
            'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
            'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
            'TwinBias_Q',E0,...
            'RingList',1:14,...
            'fixPar',fitPar,...                   % free Parameter!!
            'NonPoissonScaleFactor',1,...         % background uncertainty are enhanced
            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
            'FSDFlag','Sibille0p5eV',...          % final state distribution                        !!check ob initfile hier überschrieben wird
            'ELossFlag','KatrinT2',...            % energy loss function
            'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
            'DopplerEffectFlag','FSD',...
            'SynchrotronFlag','OFF',...
            'AngularTFFlag','OFF');
    else
        R = RunAnalysis('RunNr',10,...
            'FakeInitFile',initfile,...
            'Init_Opt',Init_Opt,...
            'RecomputeFakeRun','ON',...
            'chi2',Chi2opt,...                 % uncertainties: statistical or stat + systematic uncertainties
            'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
            'TwinBias_Q',E0,...
            'RingList',1:14,...
            'fixPar',fitPar,...                   % free Parameter!!
            'NonPoissonScaleFactor',1,...         % background uncertainty are enhanced
            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
            'FSDFlag','Sibille0p5eV',...          % final state distribution                        !!check ob initfile hier überschrieben wird
            'ELossFlag','KatrinT2',...            % energy loss function
            'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
            'DopplerEffectFlag','FSD',...
            'SynchrotronFlag','OFF',...
            'AngularTFFlag','OFF');
    end
    R.exclDataStart = R.GetexclDataStart(range);
    
    if ~isempty(Init_Opt)
        R.InitModelObj_Norm_BKG('RecomputeFlag','ON');
    end

    %% Global variables
    times = R.ModelObj.qUfrac*R.ModelObj.TimeSec;
    qU    = R.ModelObj.qU; qU    = qU-E0; % Energy axis

    % % Spectrum
    % R.ModelObj.ComputeTBDDS(); 
    % YD = R.ModelObj.TBDDS;
    % R.ModelObj.ComputeTBDIS(); 
    % YI = R.ModelObj.TBDIS; 
    % YI = YI./times;

    % Spectrum relics
    D.Fit;
    %D.ModelObj.ComputeTBDDS(); 
    YDs = D.ModelObj.TBDDS;
    %D.ModelObj.ComputeTBDIS(); 
    IS  = D.ModelObj.TBDIS; 
    YIs = IS./times;

    % Spectrum no relics
    R.Fit;
    %R.ModelObj.ComputeTBDDS(); 
    YD = R.ModelObj.TBDDS;
    %R.ModelObj.ComputeTBDIS(); 
    YI = R.ModelObj.TBDIS; 
    YI = YI./times;

    %% relic "data"

    YIsd = IS;
    % Error - stat - 
    % err  = sqrt(YIsd) ;
    % Error - stat +syst
    err  = (diag(sqrt(D.FitCMShape))) ;

    % Fluctuations (data sim)
    %YIsd = YIsd + err.*randn(length(YIsd),1);
    YIsd = YIsd./times;

    % Error bar
    err  = err./times;
    err  = err./YI;

    %% Constraining everything to qULimiteV
    qULimit = -range;
    YIsd=YIsd(qU>qULimit);
    YIs=YIs(qU>qULimit);
    YI=YI(qU>qULimit);
    sum(YI);
    err=err(qU>qULimit);
    qUc=qU(qU>qULimit);

    %% ===== PLOTTING =====
    LocalFontSize = 20;

    fig = figure('Renderer','painters');
    set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.48, 0.9]);

    plot_title = sprintf('Relic neutrino model: \\eta = %.2g',D.ModelObj.eta);
    prlG = [81 126 102]/255;
    prlB = [50 148 216]/255;
    FitStyleArg = {'o','Color','k','LineWidth',1.0,'MarkerFaceColor',rgb('Black'),'MarkerSize',4,'MarkerEdgeColor',rgb('Black')};

    % Spectra  - First Subplot

    s1=subplot(4,1,[1 2]);

    YI_N = YIs; bkg = YI_N(length(YI_N));
    YI_N = (YI_N-bkg).*(YI(1)/YI_N(1)) + bkg;

    % Plot
    %plot(qUc,YI,'DisplayName','No Sterile','color',prlB,'LineWidth',3,'LineStyle','-')
    pfit = plot(R.RunData.qU(R.exclDataStart:end)-E0,...
        R.RunData.TBDIS(R.exclDataStart:end)./...
        (R.ModelObj.qUfrac(R.exclDataStart:end)*R.ModelObj.TimeSec),...
        'DisplayName','No Relics','color',prlB,'LineWidth',3,'LineStyle','-')
    hold on

    pdata = errorbar(D.RunData.qU(D.exclDataStart:end)-E0,...
        D.RunData.TBDIS(D.exclDataStart:end)./...
        (D.ModelObj.qUfrac(D.exclDataStart:end)*D.ModelObj.TimeSec),...
        sqrt(D.RunData.TBDIS(D.exclDataStart:end))./(D.ModelObj.qUfrac(D.exclDataStart:end)*D.ModelObj.TimeSec)*50,FitStyleArg{:},'CapSize',0)

    yl1 = ylabel('Count rate (cps)');
    legend([pdata,pfit],{sprintf('Asimov data with \\eta = %.2g \n (error bars \\times 50)',D.ModelObj.eta),'\beta-decay model'},'Location','northeast','box','off');
    lgd=legend;
    lgd.FontSize = LocalFontSize-2;

    xlim([min(qUc-2) 10]);
    ylim([0.007 2*max(YI_N)]);
    PRLFormat;
    set(gca,'FontSize',LocalFontSize);
    set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
    set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
    set(gca, 'YScale', 'log');
    %yticks([1 10 100])
    if range==30
        ylim([0.007 1e2])
    else
        ylim([0.007 1e3])
    end
    %text(-41.243,31,'a)','FontSize',LocalFontSize,'FontName',get(gca,'FontName'));

    % Title
    % plot_title = sprintf('\\Deltam^{2}_{14} = %.1f eV^2 and sin^2(2\\theta_{ee}) = %.2f',sterile_mass^2,mixing_angle_2);
    % % ('Trtitum beta decay spectra comparison\nwith and without sterile neutrino\n\
    % title(plot_title)

    % Ratio    - Second Subplot

    s2=subplot(4,1,3);

    % Ratio
    RSP  = YIs./YI;
    RSPd = YIsd./YI;

    % Plot
    hr1 = plot(qUc,RSP,'color',rgb('Salmon'),'LineWidth',3,'LineStyle','-');
    hold on;
    hr2 = plot(qUc,ones(1,numel(qUc)),'color',prlB,'LineWidth',3,'LineStyle',':');
    hr3 = errorbar(qUc,RSPd,err,FitStyleArg{:},'CapSize',0);
    hold on;
    hr4 = plot(qUc,RSPd,'o','MarkerSize',2,'MarkerFaceColor',...
                    rgb('Black'),'MarkerEdgeColor',rgb('Black'),'LineWidth',3);
    pnone = plot(qUc,zeros(numel(qUc),1),'LineStyle','none');
    % Plot style
    %ylabel('Ratio \nu_4/\nu_{\beta}');
    ylabel('Ratio');
    katrinsim   = sprintf('Simulation for \\eta = %.2g',D.ModelObj.eta);
    sterilemod  = sprintf('\\beta-decay + C\\nuB model');
    hl=legend([hr2 hr1 pnone hr3],{'\beta-decay model',sterilemod,'',katrinsim},'Location','southwest','box','off');
    hl.NumColumns=2;
    hl.FontSize = LocalFontSize-2;

    xlim([min(qUc-2) 10]);

    % ylim([min([min(RSP) min(RSPd)])*0.99 max([max(RSP) max(RSPd)])*1.01])
    %ylim([min(RSP)*0.985 1.01])
    PRLFormat;
    ax2=gca;
    set(gca,'FontSize',LocalFontSize);
    set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
    set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
    %text(ax2.YLabel.Position(1)+10,1.006,'b)','FontSize',LocalFontSize,'FontName',get(gca,'FontName'));  
    ylim([0.984 1.016]);
    hl.Position(2) = 0.333;
    ax2 = gca;
    % yticks([0.9 1])

    % MTD      -  Third Subplot

    s3=subplot(4,1,4);

    bar(qUc,times(qU>qULimit)./(60*60*24),0.5,...
        'FaceColor',prlB,'EdgeColor',prlB);

    xlabel('Retarding energy - 18575 (eV)');
    ylh = ylabel('Time (d)');
    ylh.Position(1) = ax2.YLabel.Position(1);%-range-6.8;
    yl1.Position(1) = ax2.YLabel.Position(1);
    %ylim([0 230])
    %yticks([0 25 50])
    %xlim([-range-2 50])
    PRLFormat;
    set(gca,'FontSize',LocalFontSize);
    set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
    set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
    %text(ax2.YLabel.Position(1)+10,42,'c)','FontSize',get(gca,'FontSize'),'FontName',get(gca,'FontName'));      

    linkaxes([s1,s2,s3],'x');
    
    if strcmp(annotations,'ON')
       FontSize=30;
       a=annotation('doublearrow',[0.2 0.2],[0.9 0.8]);
       a.Color=[1 0 1];
       a.LineWidth=15;
       a.Head1Width=50;
       a.Head2Width=50;
       a.Head1Length=20;
       a.Head2Length=20;
       text(s1,-28.3,2.5,'N','FontSize',FontSize,'FontWeight','bold','Color',[1 0 1]);
       b=annotation('doublearrow',[0.57 0.65],[0.57 0.65]);
       b.Color=[0 0 0];
       b.LineWidth=15;
       b.Head1Width=50;
       b.Head2Width=50;
       b.Head1Length=20;
       b.Head2Length=20;
       text(s1,-12,2e-2,'m_\nu^2','FontSize',FontSize,'FontWeight','bold','Color',[0 0 0]);
       c=annotation('doublearrow',[0.65 0.75],[0.6 0.6]);
       c.Color=[1 0 0];
       c.LineWidth=15;
       c.Head1Width=50;
       c.Head2Width=50;
       c.Head1Length=20;
       c.Head2Length=20;
       text(s1,-2,0.15,'E_0^{fit}','FontSize',FontSize,'FontWeight','bold','Color',[1 0 0]);
       d=annotation('doublearrow',[0.84 0.84],[0.46 0.36]);
       d.Color=[0.9290 0.6940 0.1250];
       d.LineWidth=15;
       d.Head1Width=50;
       d.Head2Width=50;
       d.Head1Length=20;
       d.Head2Length=20;
       text(s2,8,1,'B','FontSize',FontSize,'FontWeight','bold','Color',[0.9290 0.6940 0.1250]);
       e=annotation('doublearrow',[0.665 0.665],[0.48 0.4]);
       e.Color=[81 126 102]/255;
       e.LineWidth=15;
       e.Head1Width=50;
       e.Head2Width=50;
       e.Head1Length=20;
       e.Head2Length=20;
       text(s2,-6,1.01,'\eta','FontSize',FontSize,'FontWeight','bold','Color',[81 126 102]/255);
    end

    % save
    if strcmp(saveplot,'ON')
        SaveDir = [getenv('SamakPath'),sprintf('RelicNuBkg/Plots/FinalPlots/')];
        MakeDir(SaveDir);
        SaveName='FitParams.pdf';
        export_fig(fig,[SaveDir,SaveName]);
    end
end