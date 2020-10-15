%% Settings
RunNr=1;
NetaBins=10;
etarange=10;
etafactor=1.5;
Init_Opt={'BKG_RateAllFPDSec',0.15};        % use these options by switching to RunNr 10
Recompute='OFF';

A=RelicNuDebug;
A.Chi2Scan('Recompute',Recompute,...
    'RunNr',RunNr,...
    'Init_Opt',Init_Opt,...
    'Netabins',NetaBins,...
    'etarange',etarange,...
    'etafactor',etafactor,...
    'mode','SCAN');

matFilePath = [getenv('SamakPath'),sprintf('RelicNuBkg/Chi2Scans/')];
if RunNr==1
    savename=[matFilePath,sprintf('RelicChi2Scan_Fake_NoSyst_%s_[0 %g].mat',A.Params,etafactor*10^etarange)];
elseif RunNr==10
    SaveStr='';
    for i=1:numel(Init_Opt)
        if ischar(Init_Opt{i})
            SaveStr=[SaveStr,sprintf('_%s',Init_Opt{i})];
        elseif isfloat(Init_Opt{i})
            SaveStr=[SaveStr,sprintf('_%f',Init_Opt{i})];
        end
    end
    savename=[matFilePath,sprintf('RelicChi2Scan_Fake_NoSyst_%s_[0 %g]%s.mat',A.Params,etafactor*10^etarange,SaveStr)];
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
    
    A.Chi2Scan('Recompute',Recompute,...
        'RunNr',RunNr,...
        'Init_Opt',Init_Opt,...
        'etalower',etalower,...
        'etaupper',etaupper,...
        'mode','SEARCH');
end