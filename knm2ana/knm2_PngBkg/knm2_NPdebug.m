% unblinded fit with penning track background slope
BKG_PtSlope = 3*1e-06;

NPfactor = sort(1.112+[-0.04,-0.002,-0.001,-0.0005,0,0.002,0.004,0.006,0.008,0.01,0.001,0.02]);
savedir  = [getenv('SamakPath'),'knm2ana/knm2_PngBkg/results/'];

mNuSqErrAll = zeros(numel(NPfactor),1);
mNuSqAll    = zeros(numel(NPfactor),1);
 
for i=1:numel(NPfactor)
    progressbar(i/numel(NPfactor));
    savename = sprintf('%sknm2_NPdebug_NP%.4f_BKG_PtSlope%.1fmucpsS.mat',savedir,NPfactor(i),1e6*BKG_PtSlope);
    
    if exist(savename,'file')
        d = importdata(savename);
        mNuSqErrAll(i) = d.mNuSqErr;  
         mNuSqAll(i) = d.FitResult.par(1);
        fprintf('load %s \n',savename);
    else
        SigmaSq =  0.0124+0.0025;
        range     = 40;
        RunAnaArg = {'RunList','KNM2_Prompt',...
            'DataType','Real',...
            'fixPar','mNu E0 Bkg Norm',...
            'RadiativeFlag','ON',...
            'minuitOpt','min ; minos',...
            'FSDFlag','KNM2',...
            'ELossFlag','KatrinT2A20',...
            'SysBudget',40,...
            'AnaFlag','StackPixel',...
            'chi2','chi2Stat',...
            'TwinBias_Q',18573.7,...
            'NonPoissonScaleFactor',1,...
            'FSD_Sigma',sqrt(SigmaSq),...
            'TwinBias_FSDSigma',sqrt(SigmaSq),...
            'RingMerge','Full',...
            'PullFlag',99,...
            'BKG_PtSlope',BKG_PtSlope};%99 = no pull
        
        if ~exist('A','var')
            A = MultiRunAnalysis(RunAnaArg{:});
        end
        %%
        A.exclDataStart = A.GetexclDataStart(range);
        A.NonPoissonScaleFactor = NPfactor(i);
        A.Fit;
        FitResult = A.FitResult;
        mNuSqAll(i) = A.FitResult.par(1);
        mNuSqErrAll(i) = 0.5*(A.FitResult.errPos(1)-A.FitResult.errNeg(1));
        mNuSqErr = 0.5*(A.FitResult.errPos(1)-A.FitResult.errNeg(1));
        MakeDir(savedir);
        save(savename,'FitResult','RunAnaArg','SigmaSq','BKG_PtSlope','range','mNuSqErr')
    end
end

%%
GetFigure
plot((NPfactor-1)*1e2,mNuSqErrAll,'o:','MarkerFaceColor',rgb('DodgerBlue'),'LineWidth',2)
PrettyFigureFormat;
xlabel('Non-Poisson component (%)');
ylabel(sprintf('\\sigma({\\itm}_\\nu^2) (eV^2)'));
