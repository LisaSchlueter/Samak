%-------------------------------------------------------------------------------------
% test influence of binning of neutrino mass shift induced by different RW potentials
% fit MC twins (all KNM2 runs, 3 different endpoints according to RW setting)
% replace FSD by broadened Multi-Gauss (3 peaks)
% change binning (from coarse to fine)
% Lisa Schlüter, December 2019
%-------------------------------------------------------------------------------------
% settings
E0Mode ='Runwise'; % Runwise, 3Period
RunList = 'KNM2_Prompt';
FSDFlag = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag = 'StackPixel'; % uniform FPD
chi2 = 'chi2Stat';
DataType = 'Twin';
TeStep = 0.15; % TBDDS bin size - equidistant (eV)
range = 40;
CommonArg = {'FSDFlag',FSDFlag,...
             'ELossFlag',ELossFlag,...
             'AnaFlag',AnaFlag,...
              'chi2',chi2,...
              'RunList',RunList,...
              'DataType','Twin',...
              'fixPar','mNu E0 Bkg Norm'};
%% load file
savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_TestFSDbinning_%s_%.0feVrange_%.2geV-TeBin.mat',RunList,range,TeStep)];
if exist(savename,'file')
    load(savename);
    fprintf('load %s from file \n',savename)
else
   % init Model: twins 
    switch E0Mode
        case 'Runwise'
            TwinBias_Q = 'Fit'; 
            FSDArg = {'Sigma',0.168}; %knm2 fit endpoint std
        case '3Period'
             savename3P = [savedir,'knm2_FitE0_3Periods.mat'];
            if exist(savename3P,'file')
                d = importdata(savename3P);
            else
                fprintf('file doesnt exist. run script: knm2_ApplyRWshift with E0Mode == 3Period \n')
                return
            end
             TwinBias_Q = d.TwinBias_Q;
             FSDArg = {'Sigma',0.01,'MultiGauss_RelPos',d.MultiGauss_RelPos,'MultiGauss_Weights',d.MultiGauss_Weights};
    end
    M = MultiRunAnalysis(CommonArg{:},'TwinBias_Q',TwinBias_Q);
    M.exclDataStart = M.GetexclDataStart(range);
    CommonArg = {CommonArg{:},'exclDataStart',M.exclDataStart};
    %%
    BinningFactorAll = 2:2:8;
    nBins            = numel(BinningFactorAll);
    mNuSqShift       = zeros(nBins,1);
    FSDBinningE      = cell(nBins,1);
    
    for i=1:nBins
        progressbar(i/nBins);
        M.ModelObj.SetTBDDSBinning('TeStep',TeStep);
        M.ModelObj.SetKinVariables;
        M.ModelObj.InitializeRF;
        M.ModelObj.ComputeTBDDS;
        M.ModelObj.LoadFSD(FSDArg{:},'BinningFactor',BinningFactorAll(i),'SanityPlot','OFF');
        
        M.ModelObj.ComputeTBDDS;
        M.ModelObj.ComputeTBDIS;
        M.Fit;
        mNuSqShift(i) = M.FitResult.par(1);
        FSDBinningE{i} = M.ModelObj.TTexE;
    end
    save(savename,'mNuSqShift','BinningFactorAll','nBins','TeStep','CommonArg','FSDBinningE');
    fprintf('save to file %s \n',savename)
end

%% plot
switch E0Mode
    case 'Runwise'
        f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
        p1 = plot(BinningFactorAll,mNuSqShift*1e3,'o','LineWidth',2,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
        xlabel('Bin enhancement factor');
        ylabel(sprintf('\\Delta m_\\nu^2 (meV^2)'));
        PrettyFigureFormat('FontSize',24);
        xticks(BinningFactorAll);
        xlim([BinningFactorAll(1)-0.2,BinningFactorAll(end)+0.2]);
        set(gca,'XMinorTick','off');
        savenameplot = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
      %  export_fig(f1,savenameplot);
        fprintf('save plot to %s \n',savenameplot);
end
