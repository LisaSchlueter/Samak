% Test FSD broadening with linear RW shift over time
%
RunList = 'KNM2_Prompt';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
chi2 = 'chi2Stat';
DataType = 'Twin';
freePar = 'mNu E0 Bkg Norm';
range = 40;
RunArg =     {'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2,...
    'RunList',RunList,...
    'fixPar',freePar};  % twins, all runs have same endpoint
maxRWdrift = 0.375;
savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_LinearRWShift_%s_%.0feV_maxRWdrift%.0fmeV.mat',RunList,range,maxRWdrift*1e3)];

if exist(savename,'file')
    load(savename);
else
    % Data
    D =  MultiRunAnalysis(RunArg{:},'DataType','Real'); 
    D.exclDataStart =  D.GetexclDataStart(range);
    RunArg = {RunArg{:},'exclDataStart',D.exclDataStart};
    
    %% construct distribution of RW potentials over time
    % asumption 1: overall work function change is 375 mV
    % asumption 2: work function change is linear with luminosity
    RWpotential  = maxRWdrift/D.nRuns.*(1:D.nRuns);
    TwinBias_Q = 18573.60+RWpotential;
    
    % computte twins with different endoints
    T  = MultiRunAnalysis(RunArg{:},'DataType','Twin','TwinBias_Q',TwinBias_Q); % Twins
    T.Fit; %estimate impact of averaging endpoints
    FitResultsRef = T.FitResult;
    mNuSqShift = zeros(2,1);
    mNuSqShift(1) = T.FitResult.par(1);
    
    %% manipulate FSDs: replace with rectangle distribution with width given by RW potential shift
    RecWidth = max(RWpotential)-min(RWpotential);
    FSDArg = {'Dist','Rect','Sigma',RecWidth};
    T.ModelObj.LoadFSD(FSDArg{:},'SanityPlot','OFF');
    T.ModelObj.ComputeTBDDS;
    T.ModelObj.ComputeTBDIS;
    T.Fit;
    mNuSqShift(2) = T.FitResult.par(1);
    FitResultsFSD = T.FitResult;
    save(savename,'D','T','TwinBias_Q','RWpotential','FSDArg','mNuSqShift','FitResultsRef','FitResultsFSD','RunArg');
end

%% plot plasma potential over runs
f1 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.4]);
plot(1:D.nRuns,RWpotential*1e3,'LineWidth',2.5,'Color',rgb('DodgerBlue'));
xlim([0,D.nRuns+1]);
PrettyFigureFormat('FontSize',24);
xlabel('Scan number');
ylabel(sprintf('Plasma potential change (mV)'));

plotname1 = strrep(strrep(savename,'results','plots'),'.mat','_RWpotRun.pdf');
export_fig(f1,plotname1);
%% histogram plasma potential
f2 = figure('Units','normalized','Position',[0.1,0.1,0.4,0.4]);
h1 = histogram(RWpotential*1e3,'LineWidth',0.1,'FaceColor',rgb('DodgerBlue'));
h1.BinWidth = 20;
xlim([0,D.nRuns+1]);
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('Plasma potential change (mV)'));
ylabel('Occurrence')
xlim([-10,390]);
ylim([0 21])

plotname2 = strrep(strrep(savename,'results','plots'),'.mat','RWpostHist.pdf');
export_fig(f2,plotname2);







