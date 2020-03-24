% plot neutrino mass bias as function of sigma(qU)
InitFile = @ref_FakeRun_KNM2_CD84_10hours;

range = 40;
WGTS_CD_RelStd = [0,0.05,0.1,0.2,0.5,1,2].*1e-2;
nCD = numel(WGTS_CD_RelStd);
CommonArg = {'FakeInitFile',InitFile,...
    'DataType','Fake',...
    'NonPoissonScaleFactor',1,...
    'fixPar','mNu E0 Norm Bkg'};
nRuns = 100;
savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
MakeDir(savedir);

mNuSq = zeros(nCD,1);
for i=1:nCD
    progressbar(i/nCD)
    
    savename = sprintf('%sknm2_TestStackingFakeRun_CDstd_%s_%.0feVrange_%.0fruns_%.4fRelCDStd.mat',...
        savedir, extractAfter(func2str(InitFile),'ref_FakeRun_'),...
        range,nRuns,WGTS_CD_RelStd(i));
    
    if exist(savename,'file')
        d = importdata(savename);
        mNuSq(i) = d.FitResult.par(1);
    else
        R = MultiRunAnalysis('RunList',1000+(i*nRuns:(i*nRuns+nRuns-1)),CommonArg{:},...
            'TwinBias_WGTS_CD_MolPerCm2',WGTS_CD_RelStd(i)); % in fake mode: rel. std

        R.exclDataStart = R.GetexclDataStart(range);
        R.RunData.RunName = sprintf('Stack%s_%.4fRelCDstd_%.0f_%.0f','TestCD',WGTS_CD_RelStd(i),R.RunList(1),R.RunList(end));
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
p1 =plot(qUStd,mNuSq,'o','LineWidth',2,'MarkerFaceColor',rgb('DodgerBlue'));
hold on;
x = linspace(min(qUStd),max(qUStd),100);
pref  =plot(x,-2*x.^2,'-','LineWidth',2);
pref2  =plot(x,-2*x.^2+mNuSq(2),'-.','LineWidth',2);
leg = legend([p1,pref,pref2],'361 MC runs stacked - qU randomized',...
    sprintf('{\\Delta\\itm^2} = -2\\sigma^2'),...
sprintf('{\\Delta\\itm^2} = -2\\sigma^2 + %.2f eV^2',mNuSq(2)));
leg.EdgeColor = rgb('Silver');
leg.Location='southwest';
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('\\sigma (qU)'));
ylabel(sprintf('{\\Delta\\itm^2}'));
xlim([-0.005,0.205])

saveplot = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(saveplot);