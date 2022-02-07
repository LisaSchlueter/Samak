% Choose random runlist with half of KNM1 runs
savedir = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/'];
nRunList = 1000;
FitResults = cell(nRunList,1);
RunLists   = cell(nRunList,1);
chi2 = 'chi2Stat';
NP = 1.064;

savename = [savedir,sprintf('RandomHalfRunList_Unblinded_%s_NP%2g_%.0ffits.mat',chi2,NP,nRunList)];
%%
if exist(savename,'file')
    load(savename,'FitResults','RunLists');
else
    MakeDir(savedir);
    
    RunAnaArg = { 'chi2','chi2Stat',...
    'DataType','Real',...
    'fixPar','mNu E0 Norm Bkg',... free parameter
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',NP,...
    'minuitOpt','min ; minos',...
    'FSDFlag','Sibille0p5eV',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'AngularTFFlag','OFF'};

 M = MultiRunAnalysis('RunList','KNM1_Random',RunAnaArg{:});
 
  D = copy(repmat(M,nRunList,1));
  D = reshape(D,numel(D),1);
   % progressbar('Random RunLists ');
    parfor i=1:nRunList
      %  progressbar(i/nRunList);
        D(i).RunList = D(i).GetRunList;
        D(i).nRuns         = length(D(i).RunList);
        D(i).StackRuns('CutOnSC','OFF','SCsigma',3,'CutOnFitSingleRuns','OFF');
           
        D(i).SimulateStackRuns;
        D(i).exclDataStart = D(i).GetexclDataStart(40);
        D(i).Fit;
        
      
        FitResults{i} = D(i).FitResult;
        RunLists{i} = D(i).RunList;
%        if mod(i,100)==0 % save every 100 some tmp file
%             save(strrep(savename,'.mat',sprintf('_%.0f_subsamples.mat',i)),'FitResults','RunLists');
%         end
    end
    
    save(savename,'FitResults','RunLists','M','RunAnaArg');
    fprintf('save to %s \n',savename);
end
%% histogram neutrino mass and 1 sigma
return
Both = 'OFF'; % on for blinded results
if strcmp(Both,'ON')
    savedir = [getenv('SamakPath'),'knm1ana/knm1_AlternativeRunLists/results/'];
    
    savename20 = [savedir,sprintf('RandomHalfRunList_%s_%.0f.mat',chi2,20)];
    d20 = importdata(savename20);
    mNuSq    = cell2mat(cellfun(@(x) x.par(1),d20.FitResults,'UniformOutput',0));
    mNuSqErr = cell2mat(cellfun(@(x) x.err(1),d20.FitResults,'UniformOutput',0));
    
    savename44 = [savedir,sprintf('RandomHalfRunList_%s_%.0f.mat',chi2,44)];
    d44 = importdata(savename44);
    mNuSq    = [mNuSq;cell2mat(cellfun(@(x) x.par(1),d44.FitResults,'UniformOutput',0))];
    mNuSqErr = [mNuSqErr;cell2mat(cellfun(@(x) x.err(1),d44.FitResults,'UniformOutput',0))];
     
    savename10 = [savedir,sprintf('RandomHalfRunList_%s_%.0f.mat',chi2,10)];
    d10 = importdata(savename10);
    mNuSq    = [mNuSq;cell2mat(cellfun(@(x) x.par(1),d10.FitResults,'UniformOutput',0))];
    mNuSqErr = [mNuSqErr;cell2mat(cellfun(@(x) x.err(1),d10.FitResults,'UniformOutput',0))];  
else
    mNuSq    = cell2mat(cellfun(@(x) x.par(1),FitResults,'UniformOutput',0));
    mNuSqErr = cell2mat(cellfun(@(x) x.err(1),FitResults,'UniformOutput',0));
end

PltMode = 'mNuSq';
plotdir = strrep(savedir,'results','plots');
plotname = [plotdir,sprintf('RandomHalfRunList_%s.pdf',PltMode)];

switch PltMode
    case 'mNuSq'
        x = mNuSq;
        xStr = sprintf('{\\itm}_\\nu^2 (eV^2)');
    case 'mNuSqErr'
        x = mNuSqErr;
        xStr = sprintf('\\sigma(m_\\nu^2) (eV^2)');
end

GetFigure;
h1 =histogram(mNuSq);
hold on;

h1.FaceColor = rgb('DeepSkyBlue'); h1.FaceAlpha = 1;
h1.EdgeColor = 'none';
xlabel(xStr);
ylabel('Occurence');
PrettyFigureFormat;
set(gca,'FontSize',24);
leg = legend(sprintf('%.0f samples',numel(mNuSq))); legend boxoff
%export_fig(f1,plotname);




