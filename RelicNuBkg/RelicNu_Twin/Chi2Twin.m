%% Settings
Params='KNM1';
range=40;
fitPar='mNu E0 Norm Bkg';
NetaBins=10;
etarange=11;
etafactor=5;
Recompute='OFF';

A=RelicNuDebug('Params',Params);
A.Chi2Scan_Twin('Recompute',Recompute,...
    'range',range,...
    'RunList',Params,...
    'fitPar',fitPar,...
    'Netabins',NetaBins,...
    'etarange',etarange,...
    'etafactor',etafactor,...
    'mode','SCAN');

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
savename=[matFilePath,sprintf('RelicChi2Scan_Twin_%s_[0 %g]_%s.mat',A.Params,etafactor*10^etarange,fitPar)];
                
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
        'etalower',etalower,...
        'etaupper',etaupper,...
        'mode','SEARCH');
end