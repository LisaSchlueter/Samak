function UpperLimitvsmnuSq(varargin)
    p=inputParser;
    p.addParameter('mNuBins',5,@(x)isfloat(x));
    p.addParameter('mNuUpper',4,@(x)isfloat(x));
    p.addParameter('Params','KNM1',@(x)ismember(x,{'TDR','KNM1','KNM2'}));
    p.parse(varargin{:});
    %% Settings
    mNuBins  = p.Results.mNuBins;
    mNuUpper = p.Results.mNuUpper; %eV²
    Params   = p.Results.Params;

    A = RelicNuAnalysis('Params','KNM1');

    matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
    savename = [matFilePath,sprintf('SensVsMnuSq_%s_[0_%g].mat',A.Params,mNuUpper)];

    if exist(savename,'file')
        load(savename,'ScanPoints','Sensitivities')
    else
        Sensitivities = linspace(0,mNuUpper,mNuBins);
        ScanPoints = linspace(4,mNuUpper,mNuBins);

        for i=1:numel(ScanPoints)
            A.Chi2Twin('Recompute','ON',...
                'TwinBias_mnuSq',ScanPoints(i),...
                'fitPar','mNu E0 Norm Bkg',...
                'etarange',11,...
                'etafactor',3,...
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