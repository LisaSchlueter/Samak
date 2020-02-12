thispath = 'schluete/Samak2.0';%pwd;
if contains(thispath,'Lisa/Samak2.0')         % local
    label = 'local';
elseif contains(thispath,'schluete/Samak2.0') %server
   label = 'server';
end

savedir = [getenv('SamakPath'),'knm1ana/knm1_ServerTest/results/'];
savefile = [savedir,sprintf('TestServer_FitResult_%s.mat',label)];
if exist(savefile,'file')
    load(savefile)
else
    
M = MultiRunAnalysis('RunList','KNM1','exclDataStart',17,'fixPar','5 6 7 8 9 10 11');
M.Fit
FitResult = M.FitResult;
save(savefile,'FitResult');
end