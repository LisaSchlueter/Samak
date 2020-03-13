% Compare Multiring fit with uniform fit
% MC data with radial dependence or with energy smearing
% stat. fluct
% fit once with multi-ring and once with uniform

nSamples = 1000;
Effect = 'OFF';
Random = 'ON';

if strcmp(Effect,'RadialqU')
    InitFile = @ref_FakeRun_KNM2_CD84_50days_radialqU;
elseif ismember(Effect,{'RadialmTSq','OFF'})
    InitFile = @ref_FakeRun_KNM2_CD84_50days;
end

% labeling
savedir= [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit_FakeData/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_FakeRun_Compare_MultiRing_Uniform_%s.mat',Effect)];

if strcmp(Random,'OFF')
    savename = strrep(savename,'.mat','_Asimov.mat');
elseif nSamples > 1
    savename = strrep(savename,'.mat',sprintf('_%.0fSamples.mat',nSamples));
end

if exist(savename,'file')
    load(savename)
else
    % Init model
    RunAnaArg = {'RunNr',1,...% has no meaning
        'DataType','Fake',...
        'FSDFlag','BlindingKNM2',...
        'ELossFlag','KatrinT2',...
        'exclDataStart',11,... % 11==40eV range (28 subruns)
        'chi2','chi2Stat',...
        'minuitOpt','min;migrad',...
        'NonPoissonScaleFactor',1,...
        'fitter','minuit',...
        'pullFlag',4,...
        'FakeInitFile',InitFile,....
        'fixPar','mNu E0 Bkg Norm'}; % free parameter!
    
    % Model: multi-ring fit, take qU-values from first ring
    MR = RunAnalysis(RunAnaArg{:},'AnaFlag','Ring','RingMerge','Full');
    
    if strcmp(Effect,'RadialqU')
        TBDIS_i = MR.ModelObj.TBDIS;
        qUmean = repmat(mean(MR.RunData.qU,2),[1,MR.nRings]);
        MR.SimulateRun('qU',qUmean);
    elseif strcmp(Effect,'RadialmTSq')
        MR.ModelObj.ComputeTBDDS('mTSq_bias',mTSqBias);
        MR.ModelObj.ComputeTBDIS;
        TBDIS_i= MR.ModelObj.TBDIS;
    else
        TBDIS_i= MR.ModelObj.TBDIS;
    end
    
    if strcmp(Random,'ON')
        % randomize MC data
        TBDIS = TBDIS_i+(randn([size(MR.RunData.TBDIS),nSamples]).*MR.RunData.TBDISE);
    elseif strcmp(Random,'OFF')
        TBDIS = TBDIS_i;
    end
    MR.RunData.TBDIS = TBDIS;
    
    % Uniform Model
    U = RunAnalysis(RunAnaArg{:},'AnaFlag','StackPixel');
    U.RunData.TBDIS = sum(TBDIS,2);
    
    % Fit randomized Data
    if nSamples ==1
        MR.Fit;
        U.Fit;
        Result_U = U.FitResult;
        Result_MR = MR.FitResult;
    elseif nSamples>1
        
       % Result_U  = cell(nSamples,1);
       % Result_MR = cell(nSamples,1);
        TBDIS_Sum = squeeze(sum(TBDIS,2));
        for i=1:nSamples
            progressbar(i/nSamples)
            U.RunData.TBDIS = TBDIS_Sum(:,i);
            try
            U.Fit;
            Result_U{i} = U.FitResult;
            
            
            MR.RunData.TBDIS = TBDIS(:,:,i);
            MR.Fit;
            Result_MR{i} = MR.FitResult;
            catch
            end
            if i==50
                save(savename,'Result_U','Result_MR','TBDIS','TBDIS_i','RunAnaArg');
            elseif mod(i,50)==0
                save(savename,'Result_U','Result_MR','TBDIS','TBDIS_i','RunAnaArg','-append');
            end
        end
    end
    if nSamples==1
        save(savename,'Result_U','Result_MR','TBDIS','TBDIS_i','RunAnaArg');
    end
end

fprintf('Effect: %s  - Random %s \n',Effect,Random);
if nSamples ==1
    fprintf('Uniform   mNuSq = %.3f eV \n',Result_U.par(1));
    fprintf('Multiring mNuSq = %.3f eV \n',Result_MR.par(1));
else
    fprintf('Mean Uniform   mNuSq = %.3f eV \n',mean(cell2mat(cellfun(@(x) x.par(1),Result_U,'UniformOutput',0))));
    fprintf('Mean Multiring mNuSq = %.3f eV \n',mean(cell2mat(cellfun(@(x) x.par(1),Result_MR,'UniformOutput',0))));
end
%%
par_U  = cell2mat(cellfun(@(x) x.par,Result_U,'UniformOutput',0))';
par_MR = cell2mat(cellfun(@(x) x.par,Result_MR,'UniformOutput',0))';

% plot histogram
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
hU = histogram(par_U(1,:));
hU.FaceAlpha = 0.5;
hU.FaceColor = rgb('DodgerBlue');
hU.EdgeColor = rgb('DodgerBlue');
hold on;
hMR = histogram(par_MR(1,:));
hMR.FaceColor = rgb('ForestGreen');
hMR.FaceAlpha = 0.5;
hMR.EdgeColor = rgb('ForestGreen');%'none';
thiYlim = ylim;
%plot(mean(par_U(1,:)).*ones(2,1),[0,max(thiYlim)],'Color',hU.FaceColor,'LineWidth',2);
%plot(mean(par_MR(1,:)).*ones(2,1),[0,max(thiYlim)],'Color',hMR.FaceColor,'LineWidth',2);
leg = legend([hU,hMR],'Uniform fit','Multi-Ring fit');
legend boxoff;
leg.Location = 'northwest';
PrettyFigureFormat('FontSize',24);
xstr = sprintf('{\\itm}_\\nu^2 (eV^2)');
xlabel(xstr);
ylabel('Events');
savedirplot = strrep(savedir,'results','plots');
MakeDir(savedirplot);
savenameplot1 = strrep(strrep(savename,'results','plots'),'.mat','_Hist.pdf');
export_fig(f1,savenameplot1,'-painters');
%% histogram difference
f2 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
hDiff = histogram(par_U(1,:)-par_MR(1,:));
hDiff.FaceAlpha = 0.8;
hDiff.FaceColor = rgb('DodgerBlue');
hDiff.EdgeColor = rgb('SteelBlue');
thiYlim = ylim;
MeanDiff = mean(par_U(1,:)-par_MR(1,:));
hold on;
pMean = plot(MeanDiff.*ones(2,1),[0,max(thiYlim)],'Color','k','LineWidth',2);
%plot(mean(par_MR(1,:)).*ones(2,1),[0,max(thiYlim)],'Color',hMR.FaceColor,'LineWidth',2);
leg = legend([hDiff,pMean],sprintf('\\DeltaFit result (U - MR)'),sprintf('Mean = %.3f eV^2',MeanDiff));
legend boxoff;
leg.Location = 'northwest';
PrettyFigureFormat('FontSize',24);
xstr = sprintf('\\Delta{\\itm}_\\nu^2 (eV^2)');
xlabel(xstr);
ylabel('Events');

savedirplot = strrep(savedir,'results','plots');
MakeDir(savedirplot);
savenameplot2 = strrep(strrep(savename,'results','plots'),'.mat','_HistDiff.pdf');
export_fig(f2,savenameplot2,'-painters');
