% extend grid for those randMC that yielded best fit on grid boundary

% grid search on randomized twins
% ksn2 calculate chi2 grid search
randMC = 1:1e3;
NrandMC = 1e3;
Twin_sin2T4 = 0;
Twin_mNu4Sq = 0;
chi2 = 'chi2CMShape';
DataType = 'Twin';
freePar = 'E0 Norm Bkg';
nGridSteps = 25;
range = 40;
Mode ='Interp';% 'Compute';
%% get summary file
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_%.0fsamples.mat',savedir,NrandMC);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,NrandMC);
end

if exist(savefile,'file') 
    fprintf('load file %s \n',savefile);
    d = importdata(savefile);
else
     fprintf('file missing %s \n',savefile);
    return
end

% find fits on grid boundary: sin2t4==0.5
Idx= find(d.sin2T4_bf==0.5);

%% configure RunAnalysis object
if strcmp(chi2,'chi2Stat')
    NonPoissonScaleFactor = 1;
elseif  strcmp(chi2,'chi2CMShape')
    NonPoissonScaleFactor = 1.112;
end
RunAnaArg = {'RunList','KNM2_Prompt',...
    'DataType',DataType,...
    'fixPar',freePar,...%free par
    'SysBudget',40,...
    'fitter','minuit',...
    'minuitOpt','min;migrad',...
    'RadiativeFlag','ON',...
    'FSDFlag','KNM2_0p5eV',...
    'ELossFlag','KatrinT2A20',...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'NonPoissonScaleFactor',NonPoissonScaleFactor,...
    'FSD_Sigma',sqrt(0.0124+0.0025),...
    'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
    'TwinBias_Q',18573.7,...
    'PullFlag',99,...;%99 = no pull
    'BKG_PtSlope',3*1e-06,...
    'TwinBias_BKG_PtSlope',3*1e-06,...
    'DopplerEffectFlag','FSD'};
A = MultiRunAnalysis(RunAnaArg{:});
%% configure Sterile analysis object
SterileArg = {... % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    'nGridSteps',nGridSteps,...
    'SmartGrid','OFF',...
    'RecomputeFlag','OFF',...
    'SysEffect','all',...
    'range',range,...
    'RandMC','OFF',...
    'Twin_mNu4Sq',Twin_mNu4Sq,...
    'Twin_sin2T4',Twin_sin2T4};
%%
switch Mode
    case 'Compute'
        for i=1:numel(Idx)
            A = MultiRunAnalysis(RunAnaArg{:});
            S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
            S.RandMC= randMC(Idx(i));
            %% load file to get tbdis_mc
            df = importdata(S.GridFilename('mNu4SqTestGrid',2,'ExtmNu4Sq','ON'));
            S.RandMC_TBDIS = df.TBDIS_mc;
            S.GridSearch('mNu4SqTestGrid',6,'ExtmNu4Sq','ON');
        end
    case 'Interp'
        A = MultiRunAnalysis(RunAnaArg{:});
        S = SterileAnalysis('RunAnaObj',A,SterileArg{:});
        %%
        i = 1;
     %   for i=1:numel(Idx)
            % load file to get tbdis_mc
            S.RandMC_TBDIS = [];
            S.RandMC= randMC(Idx(i));
            filename1 = S.GridFilename('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
            df = importdata(filename1);
           
            % test
            S.LoadGridFile('mNu4SqTestGrid',2,'ExtmNu4Sq','ON');
            S.Interp1Grid('Minm4Sq',0.1,'Maxm4Sq',5);
            mNu4Sq = S.mNu4Sq;
            sin2T4 = S.sin2T4;
            chi2 = S.chi2;
            %end
            
            S.RandMC_TBDIS = df.TBDIS_mc;
            S.LoadGridFile('mNu4SqTestGrid',6,'ExtmNu4Sq','ON');
            S.Interp1Grid('Minm4Sq',0.1,'Maxm4Sq',5);
            mNu4Sq_2 = S.mNu4Sq;
            sin2T4_2 = S.sin2T4;
            chi2_2 = S.chi2;
            
            mNu4Sq_tot = [mNu4Sq;mNu4Sq_2];
            sin2T4_tot = [sin2T4;sin2T4_2];
            chi2_tot = [chi2;chi2_2];
            
            GetFigure;
            s1 =surf(sin2T4_tot,mNu4Sq_tot,chi2_tot,'EdgeColor','interp','FaceColor','interp');
            view(2)
            set(gca,'XScale','lin')
            set(gca,'YScale','log')
            xlim([1e-03 1])
            ylim([0.1 5])
            c = colorbar;
          
            %% cross check chi2 of 1 grid point
            A.RunData.TBDIS = df.TBDIS_mc;
            A.ModelObj.SetFitBiasSterile(1,0.5)
            A.Fit;
            %% end
%             S.FindBestFit
%             BestFit_Ext = struct('chi2_bf',S.chi2_bf,'mNu4Sq_bf',S.mNu4Sq_bf,...
%                 'sin2T4_bf',S.sin2T4_bf,'mNuSq_bf',S.mNuSq_bf,'E0_bf',S.E0_bf);
%          %   save(filename1,'BestFit_Ext','-append');
%         %    fprintf('save extended best fit to %s \n',filename1);
%             
%       %  end
end
