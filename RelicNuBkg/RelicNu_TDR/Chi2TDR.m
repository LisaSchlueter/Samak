%% Settings
Params='KNM1';
range=50;
RunNr=10;
NetaBins=10;
etarange=11;
etafactor=1;
fitPar='mNu E0 Norm Bkg';
Init_Opt={'BKG_RateAllFPDSec',0.4,'mNuSq_i',1,'TimeSec',3*365.242*24*3600};        % use these options by switching to RunNr 10
Syst='OFF';
Recompute='OFF';

A=RelicNuDebug('Params',Params);
A.Chi2Scan_Fake('Recompute',Recompute,...
    'range',range,...
    'RunNr',RunNr,...
    'fitPar',fitPar,...
    'Syst',Syst,...
    'Init_Opt',Init_Opt,...
    'Netabins',NetaBins,...
    'etarange',etarange,...
    'etafactor',etafactor,...
    'mode','SCAN');

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
if RunNr==1
    savename=[matFilePath,sprintf('RelicChi2Scan_Fake_Syst%s_range%g_%s_[0 %g]_%s.mat',Syst,range,Params,etafactor*10^etarange,fitPar)];
elseif RunNr==10
    SaveStr='';
    for i=1:numel(Init_Opt)
        if ischar(Init_Opt{i})
            SaveStr=[SaveStr,sprintf('_%s',Init_Opt{i})];
        elseif isfloat(Init_Opt{i})
            SaveStr=[SaveStr,sprintf('_%f',Init_Opt{i})];
        end
    end
    savename=[matFilePath,sprintf('RelicChi2Scan_Fake_Syst%s_range%g_%s_[0 %g]%s_%s.mat',Syst,range,Params,etafactor*10^etarange,SaveStr,fitPar)];
end
                
load(savename);

if Chi2(end)<2.71
    sprintf('Increase etafactor or etarange')
else
    etavalues=(0:(Netabins-1))*((etafactor*10^(etarange))/(Netabins-1));
    m=find(Chi2<2.71);
    n=find(Chi2>2.71);
    etalower=etavalues(m(end));
    etaupper=etavalues(n(1));
    
    A.Chi2Scan_Fake('Recompute',Recompute,...
        'range',range,...
        'RunNr',RunNr,...
        'fitPar',fitPar,...
        'Syst',Syst,...
        'Init_Opt',Init_Opt,...
        'etalower',etalower,...
        'etaupper',etaupper,...
        'mode','SEARCH');
end