function Chi2Twin(varargin)
    p=inputParser;
    p.addParameter('Params','KNM1',@(x)ismember(x,{'TDR','KNM1','KNM2'}));
    p.addParameter('range',40,@(x)isfloat(x));
    p.addParameter('fitPar','mNu E0 Norm Bkg',@(x)ischar(x));
    p.addParameter('Syst','ON',@(x)ismember(x,{'ON','OFF'}));
    p.addParameter('SystBudget',24,@(x)isfloat(x));
    p.addParameter('TwinBias_mnuSq',1,@(x)isfloat(x));
    p.addParameter('NetaBins',10,@(x)isfloat(x));
    p.addParameter('etarange',11,@(x)isfloat(x));
    p.addParameter('etafactor',5,@(x)isfloat(x));
    p.addParameter('Recompute','OFF',@(x)ismember(x,{'ON','OFF'}));
    p.parse(varargin{:});
    
    Params         = p.Results.Params;
    range          = p.Results.range;
    fitPar         = p.Results.fitPar;
    Syst           = p.Results.Syst;
    SystBudget     = p.Results.SystBudget;
    TwinBias_mnuSq = p.Results.TwinBias_mnuSq;
    NetaBins       = p.Results.NetaBins;
    etarange       = p.Results.etarange;
    etafactor      = p.Results.etafactor;
    Recompute      = p.Results.Recompute;
    

    A=RelicNuDebug('Params',Params);
    A.Chi2Scan_Twin('Recompute',Recompute,...
        'range',range,...
        'RunList',Params,...
        'fitPar',fitPar,...
        'Syst',Syst,...
        'SystBudget',SystBudget,...
        'TwinBias_mnuSq',TwinBias_mnuSq,...
        'Netabins',NetaBins,...
        'etarange',etarange,...
        'etafactor',etafactor,...
        'mode','SCAN');

    matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
    savename=[matFilePath,sprintf('RelicChi2Scan_Twin_BiasmnuSq%g_Syst%s_range%g_%s_[0 %g]_%s.mat',TwinBias_mnuSq,Syst,range,Params,etafactor*10^etarange,fitPar)];

    load(savename);

    if Chi2(end)<2.71
        sprintf('Increase etafactor or etarange')
    else
        etavalues=(0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1));
        m=find(Chi2<2.71);
        n=find(Chi2>2.71);
        etalower=etavalues(m(end));
        etaupper=etavalues(n(1));

        A.Chi2Scan_Twin('Recompute',Recompute,...
            'range',range,...
            'RunList',Params,...
            'fitPar',fitPar,...
            'Syst',Syst,...
            'SystBudget',SystBudget,...
            'TwinBias_mnuSq',TwinBias_mnuSq,...
            'etalower',etalower,...
            'etaupper',etaupper,...
            'mode','SEARCH');
    end
end