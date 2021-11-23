function SensitivityVsmnuSq(varargin)
    p=inputParser;
    p.addParameter('mNuBins',5,@(x)isfloat(x));
    p.addParameter('mNuUpper',4,@(x)isfloat(x));
    p.addParameter('Params','TDR',@(x)ismember(x,{'TDR','KNM1','KNM2'}));
    p.addParameter('T','OFF',@(x)ismember(x,{'ON','OFF'}));             %Atomic Tritium
    p.addParameter('TDRbkg','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('MTD','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('DeltaChi2',2.71,@(x)isfloat(x));
    p.parse(varargin{:});
    %% Settings
    mNuBins  = p.Results.mNuBins;
    mNuUpper = p.Results.mNuUpper; %eV²
    Params   = p.Results.Params;
    T        = p.Results.T;
    TDRbkg   = p.Results.TDRbkg;
    MTD      = p.Results.MTD;
    DeltaChi2= p.Results.DeltaChi2;

    A = RelicNuAnalysis('Params',Params);

    matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/')];
    savename = [matFilePath,sprintf('SensVsmnuSq_%s_[0_%g]',A.Params,mNuUpper)];
    if strcmp(T,'ON')
        savename = [savename,'_T'];
    end
    if strcmp(TDRbkg,'ON')
        savename = [savename,'_TDRbkg'];
    end
    if strcmp(MTD,'ON')
        savename = [savename,'_MTD'];
    end
    if strcmp(T,'OFF') && strcmp(TDRbkg,'OFF') && strcmp(MTD,'OFF')
        savename = [savename,'_Syst'];
    end
    savename = [savename,'.mat'];

    if exist(savename,'file')
        load(savename,'ScanPoints','Sensitivities')
    else
        Sensitivities = zeros(1,mNuBins);
        ScanPoints = linspace(0,mNuUpper,mNuBins);
        if strcmp(Params,'TDR')
            etafactor=2;
            etarange=10;
        else
            etafactor=3;
            etarange=11;
        end

        for i=1:numel(ScanPoints)
            Init_Opt={'mnuSq_i',ScanPoints(i)};
            if strcmp(T,'ON')
                Init_Opt = [Init_Opt,{'TTFSD','OFF','HTFSD','OFF','DTFSD','OFF'}];
            end
            if strcmp(TDRbkg,'OFF')
                Init_Opt = [Init_Opt,{'BKG_RateAllFPDSec',0.13}];
            end
            if strcmp(T,'OFF') && strcmp(TDRbkg,'OFF') && strcmp(MTD,'OFF')
                Syst='ON';
            else
                Syst='OFF';
            end
            A.Chi2Fake('Recompute','ON',...
                'RunNr',10,...
                'range',40,...
                'etafactor',etafactor,...
                'etarange',etarange,...
                'fitPar','mNu E0 Norm Bkg',...
                'Syst',Syst,...
                'Init_Opt',Init_Opt,... %,'TTFSD','OFF','HTFSD','OFF','DTFSD','OFF'},...
                'NetaBins',2,...
                'Plot','OFF',...
                'DeltaChi2',DeltaChi2);
            Sensitivities(i) = A.etaSensitivity;
        end

        save(savename,'ScanPoints','Sensitivities');
    end
    fig1=figure(1);
plot(ScanPoints,Sensitivities,'LineWidth',2);
xlabel('m_\nu^2','FontSize',12);
ylabel('\eta','FontSize',12);
PrettyFigureFormat;
end

