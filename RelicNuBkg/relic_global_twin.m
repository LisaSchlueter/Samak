function relic_global_twin(varargin)

    %% Settings
    p=inputParser;
    p.addParameter('eta',2.8e9,@(x)isfloat(x));
    p.addParameter('RunList','KNM1',@(x)ismember(x,{'KNM1','KNM2'}));
    p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
    p.addParameter('E0',18573.73,@(x)isfloat(x));                                         % Endpoint in eV
    p.addParameter('range',40,@(x)isfloat(x));
    p.parse(varargin{:});
    eta      =p.Results.eta;
    RunList  =p.Results.RunList;
    fitPar   =p.Results.fitPar;
    E0       =p.Results.E0;
    range    =p.Results.range;
    
    switch RunList
        case 'KNM2'
            initfile=@ref_RelicNuBkg_KNM2;
        case 'KNM1'
            initfile=@ref_RelicNuBkg_KNM1;
    end

    %% Data
    D = RunAnalysis('RunNr',100,...
        'RecomputeFakeRun','ON',...
        'Init_Opt',{'eta_i',eta},...
        'FakeInitFile',initfile,...
        'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Fake',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'TwinBias_Q',E0,...
        'fixPar',fitPar,...                   % free Parameter!!
        'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...           % final state distribution                        !!check ob initfile hier überschrieben wird
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'DopplerEffectFlag','FSD',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF');
    
    %D.InitModelObj_Norm_BKG('RecomputeFlag','ON');

    D.exclDataStart = D.GetexclDataStart(range);
    %% No relics

    if strcmp(RunList,'KNM1')
        R = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
            'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
            'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
            'fixPar',fitPar,...                   % free Parameter!!
            'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
            'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
            'FSDFlag','SibilleFull',...           % final state distribution
            'ELossFlag','KatrinT2',...            % energy loss function
            'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
            'DopplerEffectFlag','FSD',...
            'Twin_SameCDFlag','OFF',...
            'Twin_SameIsotopFlag','OFF',...
            'SynchrotronFlag','ON',...
            'AngularTFFlag','OFF',...
            'TwinBias_Q',18573.73);
    elseif strcmp(RunList,'KNM2')
        R = MultiRunAnalysis('RunList',RunList,... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
            'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
            'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
            'fixPar',fitPar,...                   % free Parameter!!
            'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
            'NonPoissonScaleFactor',1.064,...     % background uncertainty are enhanced
            'minuitOpt','min ; minos',...         % technical fitting options (minuit)
            'FSDFlag','SibilleFull',...           % final state distribution
            'ELossFlag','KatrinT2',...            % energy loss function
            'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
            'DopplerEffectFlag','FSD',...
            'Twin_SameCDFlag','OFF',...
            'Twin_SameIsotopFlag','OFF',...
            'SynchrotronFlag','ON',...
            'AngularTFFlag','OFF',...
            'TwinBias_Q',18573.73);
    end
    R.exclDataStart = R.GetexclDataStart(range);

    %% Global variables
    times = R.ModelObj.qUfrac*R.ModelObj.TimeSec;
    qU    = R.ModelObj.qU; qU    = qU-E0; % Energy axis

    % Spectrum no relics
    R.Fit;
    YD = R.ModelObj.TBDDS;
    YI = R.ModelObj.TBDIS; 
    YI = YI./times;
    
    % Spectrum relics
    D.ModelObj.NormFactorTBDDS=(R.FitResult.par(3+R.ModelObj.nPixels:3+2*R.ModelObj.nPixels-1) + 1).*R.ModelObj.NormFactorTBDDS;
    D.ModelObj.ComputeTBDDS(); 
    YDs = D.ModelObj.TBDDS;
    D.ModelObj.ComputeTBDIS(); 
    IS  = D.ModelObj.TBDIS;
    YIs = IS./times;

    %% relic "data"

    YIsd = IS;
    % Error - stat - 
    % err  = sqrt(YIsd) ;
    % Error - stat +syst
    err  = (diag(sqrt(R.FitCMShape))) ;

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
    set(fig, 'Units', 'normalized', 'Position', [0.001, 0.001,0.45, 0.8]);

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
        (D.ModelObj.TBDIS(D.exclDataStart:end))./...
        (D.ModelObj.qUfrac(D.exclDataStart:end)*D.ModelObj.TimeSec),...
        sqrt(D.RunData.TBDIS(D.exclDataStart:end))./(D.ModelObj.qUfrac(D.exclDataStart:end)*D.ModelObj.TimeSec)*50,FitStyleArg{:},'CapSize',0)

    yl1 = ylabel('Count rate (cps)');
    legend([pdata,pfit],{sprintf('Asimov data with \\eta = %.2g (error bars \\times 50)',D.ModelObj.eta),'\beta-decay model'},'Location','northeast','box','off');
    lgd=legend;
    lgd.FontSize = LocalFontSize-2;

    xlim([min(qUc-5) 10]);
    ylim([0.007 2*max(YI_N)]);
    PRLFormat;
    set(gca,'FontSize',LocalFontSize);
    set(get(gca,'XLabel'),'FontSize',LocalFontSize+4);
    set(get(gca,'YLabel'),'FontSize',LocalFontSize+4);
    set(gca, 'YScale', 'log');
    %yticks([1 10 100])
    ylim([0.007 1e3])
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

    xlim([min(qUc-5) 10]);

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

    % save
    export_fig(fig,'./plots/relic_spectrum.pdf');
end