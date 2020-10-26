%% Settings
Params='KNM1';
range=40;
fitPar='mNu E0 Norm Bkg';
Syst='OFF';
TwinBias_mnuSq=0;
NetaBins=10;
etarange=11;
etafactor=5;
Recompute='OFF';

A=RelicNuDebug('Params',Params);
A.Chi2Scan_Twin('Recompute',Recompute,...
    'range',range,...
    'RunList',Params,...
    'fitPar',fitPar,...
    'Syst',Syst,...
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
        'TwinBias_mnuSq',TwinBias_mnuSq,...
        'etalower',etalower,...
        'etaupper',etaupper,...
        'mode','SEARCH');
end