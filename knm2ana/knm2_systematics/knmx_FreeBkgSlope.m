% background slope as a free fit parameter for final KATRIN settings
range = 40;
InitFile = @ref_FakeRun_FinalKATRIN_CD100_1000days;
NPfactor = 1;%.112;
savedir = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2xFreeBkgSlope.mat',savedir);

if exist(savename,'file')
    load(savename)
else
%% init model
F = RunAnalysis('RunNr',1,'FakeInitFile',InitFile,...
    'DataType','Fake',...
    'fixPar','mNu E0 Bkg Norm',...
    'NonPoissonScaleFactor',NPfactor,...
    'PixList',1:124,...
    'FSDFlag','Sibille0p5eV',...
    'ElossFlag','KatrinT2A20',...
    'ROIFlag','Default');

F.exclDataStart = F.GetexclDataStart(range);

%% stat. only - fixed background slope
F.Fit;
FitResults_BsFix = F.FitResult;

%% stat. only - free background slope
F.fixPar = 'mNu E0 Bkg Norm BkgSlope';
F.InitFitPar;
F.Fit;
FitResults_BsFree = F.FitResult;
save(savename,'FitResults_BsFix','FitResults_BsFree','range','InitFile','NPfactor');
end

%% result
mNuSqErrFree = 0.5*(-FitResults_BsFree.errNeg(1)+FitResults_BsFree.errPos(1));
mNuSqErrFix  = 0.5*(-FitResults_BsFix.errNeg(1)+FitResults_BsFix.errPos(1));

fprintf('----------------------------------------------\n');
fprintf('nu-mass sq. syst. sensitivity: free slope = %.3f eV^2 \n',sqrt(mNuSqErrFree^2-mNuSqErrFix^2));
fprintf('background slope sensitivity  = %.2f mcps / keV\n',FitResults_BsFree.err(12)*1e6);

