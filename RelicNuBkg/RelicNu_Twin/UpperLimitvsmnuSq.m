function UpperLimitvsmnuSq(varargin)
    p=inputParser;
    p.addParameter('mNuBins',5,@(x)isfloat(x));
    p.addParameter('mNuUpper',4,@(x)isfloat(x));
    p.addParameter('Params','KNM1',@(x)ismember(x,{'TDR','KNM1','KNM2_Prompt'}));
    p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    %% Settings
    mNuBins  = p.Results.mNuBins;
    mNuUpper = p.Results.mNuUpper; %eV²
    Params   = p.Results.Params;
    Recompute= p.Results.Recompute;

    A = RelicNuAnalysis('Params',Params);
    if strcmp(Params,'KNM1')
        SB=24;
        etafactor=3;
        etarange=11;
    elseif strcmp(Params,'KNM2_Prompt')
        SB=40;
        etafactor=1.3;
        etarange=11;
    end

    matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
    savename = [matFilePath,sprintf('SensVsMnuSq_%s_[0_%g].mat',A.Params,mNuUpper)];

    if exist(savename,'file')
        load(savename,'ScanPoints','Sensitivities')
    else
        Sensitivities = linspace(0,mNuUpper,mNuBins);
        ScanPoints = linspace(0,mNuUpper,mNuBins);

        for i=1:numel(ScanPoints)
            A.Chi2Twin('Recompute',Recompute,...
                'TwinBias_mnuSq',ScanPoints(i),...
                'fitPar','mNu E0 Norm Bkg',...
                'SystBudget',SB,...
                'etarange',etarange,...
                'etafactor',etafactor,...
                'NetaBins',2,...
                'DeltaChi2',2.71,...
                'Plot','OFF');
            Sensitivities(i) = A.etaSensitivity;
        end

        save(savename,'ScanPoints','Sensitivities');
    end

    fig1=figure(1);
    plot(ScanPoints,Sensitivities,'LineWidth',2);
    xlabel('m_{\nu}^{2} (eV^{2})','FontSize',12);
    ylabel('\eta','FontSize',12);
    PrettyFigureFormat;
end