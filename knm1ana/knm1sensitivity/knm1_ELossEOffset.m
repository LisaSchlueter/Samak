RunList = 'KNM1';
DataType = 'Twin';
exclDataStart = 14;
chi2 = 'chi2Stat';
is_EOffset = (-80:10:80)*1e-03;
RecomputeFlag = 'ON';

nE = numel(is_EOffset);
M = MultiRunAnalysis('RunList',RunList,'DataType',DataType,'exclDataStart',exclDataStart,'fixPar','5 6 7 8 9 10 11',...
        'SysBudget',2,'chi2',chi2);
    
savedir = [getenv('SamakPath'),'knm1ana/knm1_ElossOffset/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('ElossOffset_%s_%s_%s_%.0f_minE%.3fmeV_maxE%.3fmeV.mat',RunList,DataType,chi2,exclDataStart,...
    min(is_EOffset)*1e3,1e3*max(is_EOffset))];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename,mNuSq,mNuSqErr);
else
  
    mNuSq    = zeros(nE,1);
    mNuSqErr = zeros(nE,1);
    chi2min  = zeros(nE,1);
    
    for i=1:nE
        M.ModelObj.is_EOffset = is_EOffset(i);
        M.ModelObj.InitializeRF;
        M.Fit;
        mNuSq(i)     = M.FitResult.par(1);
        mNuSqErr(i)  = M.FitResult.err(1);
        chi2min(i)   = M.FitResult.chi2min;
    end
    save(savename,'mNuSq','mNuSqErr','chi2min','is_EOffset');
end
%%
plot(is_EOffset*1e3,mNuSq,'x--','LineWidth',3);
legend(sprintf('%.0f eV range',M.RunData.qU(exclDataStart)-18573.7));
legend boxoff; 
PrettyFigureFormat;