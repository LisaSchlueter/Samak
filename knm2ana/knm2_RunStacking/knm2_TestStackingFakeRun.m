
InitFile = @ref_FakeRun_KNM2_CD84_8hours;
range = 40;
TestOpt = {'AllSame','qUfrac','qU'};%,'qUqUfrac'};
CommonArg = { 'FakeInitFile',InitFile,...
    'DataType','Fake',...
    'NonPoissonScaleFactor',1,...
    'fixPar','mNu E0 Norm Bkg'};
nRuns = 100;
savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
MakeDir(savedir);

mNuSq = zeros(numel(TestOpt),1);

for i=1:numel(TestOpt)
    progressbar(i/numel(TestOpt))
    
    savename = sprintf('%sknm2_TestStackingFakeRun_%s_%.0feVrange_%.0fruns_%s.mat',...
        savedir, extractAfter(func2str(InitFile),'ref_FakeRun_'),...
        range,nRuns,TestOpt{i});
    
    if exist(savename,'file')
        d = importdata(savename);
        mNuSq(i) = d.FitResult.par(1);
    else
        
        if strcmp(TestOpt{i},'AllSame')
            R = MultiRunAnalysis('RunList',1000+(i*nRuns:(i*nRuns+nRuns-1)),CommonArg{:});% default
        elseif strcmp(TestOpt{i},'qUfrac')
            R = MultiRunAnalysis('RunList',1000+(i*nRuns:(i*nRuns+nRuns-1)),CommonArg{:},...
                'TwinBias_qUfrac',0.2); % in fake mode: relative std
        elseif strcmp(TestOpt{i},'qU')
            R = MultiRunAnalysis('RunList',1000+(i*nRuns:(i*nRuns+nRuns-1)),CommonArg{:},...
                'TwinBias_qU',0.01); % in fake mode: abs std
        elseif strcmp(TestOpt{i},'qUqUfrac')
            R = MultiRunAnalysis('RunList',1000+(i*nRuns:(i*nRuns+nRuns-1)),CommonArg{:},...
                'TwinBias_qU',0.01,... % in fake mode: abs std
                'TwinBias_qUfrac',0.2);
        end
        
        R.exclDataStart = R.GetexclDataStart(range);
        R.RunData.RunName = sprintf('Stack%s_%.0f_%.0f',TestOpt{i},R.RunList(1),R.RunList(end));
        R.Fit;
        mNuSq(i) = R.FitResult.par(1);
        
        %save
        FitResult = R.FitResult;
        TBDIS_d = R.RunData.TBDIS;
        TBDIS_m = R.ModelObj.TBDIS;
        qU = R.ModelObj.qU;
        SingleRunData = R.SingleRunData;
        RunData       = R.RunData;
        save(savename,'FitResult','qU','TBDIS_d','TBDIS_m','SingleRunData','RunData');
    end  
end
%% plot
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
p1 =plot(1:numel(TestOpt),mNuSq,':o','LineWidth',2);
hold on;
leg = legend('100 MC runs stacked','Location','northwest');
leg.EdgeColor = rgb('Silver');
PrettyFigureFormat('FontSize',24);
xticks(1:numel(TestOpt));
xticklabels({'all same','rand qUfrac',sprintf('rand qU \\sigma = 10 meV')});
ylabel(sprintf('{\\Delta\\itm^2}'));
xlim([0.8 3.2]);

saveplot = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(saveplot);