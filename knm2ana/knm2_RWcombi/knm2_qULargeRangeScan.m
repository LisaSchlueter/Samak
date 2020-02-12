% idea: infer RW shift from data, which is not used in fit
% test uniform FPD (and ringwise later)
% fit E0 as function of upper fit boundary [-90, -40],  [-90, -30] , etc.
% goal: check if variation small and that the large range fit is realiable
% result:
RecomputeFlag = 'OFF';
RunList   = 'KNM2_RW3';%'KNM2_Prompt';
FSDFlag   = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag   = 'StackPixel'; % uniform FPD
chi2      = 'chi2Stat';
freePar   = 'E0 Norm'; %free fit parameter
DataType  = 'Real';

RunArg =     {'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2,...
    'RunList',RunList,...
    'fixPar',freePar}; % all twins have same endpoint


savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_qULargeRangeScan_%s_%s_freePar-%s.mat',DataType,RunList,strrep(freePar,' ',''))];

%% load
if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    %% RunAnalysis Object: Real
    T =  MultiRunAnalysis(RunArg{:});
    
    % define analysis energy range
    T.exclDataStart = 1;                               % start at first bin: deepest into spectrum
    T.exclDataStop =  T.GetexclDataStart(40)-1;        % stop at 1 before 40eV
    RunArg = {RunArg{:},'exclDataStart',T.exclDataStart};
    
    %%
    BkgMean = mean(T.RunData.TBDIS(T.RunData.qU>18575)./(T.RunData.qUfrac(T.RunData.qU>18575).*T.RunData.TimeSec));
    T.ModelObj.BKG_RateSec_i = BkgMean;
    T.ModelObj.SetFitBias(0);
    
    ranges = [40,30,20 10];
    exclDataStop = zeros(numel(ranges),1);
    par = zeros(T.nPar,numel(ranges));
    err = zeros(T.nPar,numel(ranges));
    pval = zeros(numel(ranges),1);
    
    
    for i=1:numel(ranges)
        progressbar(i/numel(ranges));
        T.exclDataStop =  T.GetexclDataStart(ranges(i))-1;
        
        T.Fit;
        par(:,i) = T.FitResult.par;
        err(:,i) = T.FitResult.err;
        pval(i)  = 1-chi2cdf(T.FitResult.chi2min,T.FitResult.dof);
        exclDataStop(i) = T.exclDataStop;
    end
    
    upperfitrange = T.RunData.qU(exclDataStop)-T.ModelObj.Q_i;
    RunData = T.RunData;
    Q_i = T.ModelObj.Q_i;
    save(savename,'RunData','Q_i','par','err','ranges','upperfitrange','pval','exclDataStop','BkgMean','RunArg','T');
end
%% plot
f1 = figure('Units','normalized','Position',[0.1,0.1,0.45,0.45]);
y = (par(2,:)-mean(par(2,:))).*1e3;
yErr = err(2,:).*1e3;
[l,a] = boundedline(linspace(RunData.qU(11)-Q_i,0,10),zeros(10,1),100.*ones(10,1));
a.FaceColor = rgb('LightGray'); a.FaceAlpha = 0.7;
l.LineStyle = 'none';
hold on;
l1 = plot(linspace(min(upperfitrange)-2,max(upperfitrange)+2,10),zeros(10,1),'--','Color',rgb('Silver'),'LineWidth',2);
hold on
e1 = errorbar(upperfitrange,y,yErr,'-o','CapSize',0,'LineWidth',l1.LineWidth,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'));
hold off;
PrettyFigureFormat('FontSize',24);
xlabel('Upper fit boundary below E_0 (eV)');
ylabel(sprintf('{\\itE}_0 - \\langle{\\itE}_0\\rangle (meV)'))
xlim([min(upperfitrange)-2,max(upperfitrange)+2])
ylim([min(y-yErr)*1.2,max(y+yErr).*1.2])
leg = legend([e1,a],sprintf(' Free fit parameter: %s \n Lower fit boundary: -90 eV',freePar),sprintf(' Data used for m_\\nu fit'));%,sprintf('\\langle{\\itE}_0\\rangle'));
leg.Location='southeast';
legend boxon;
leg.EdgeColor = 'none';
leg.Color = rgb('White');

% save plot
plotdir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/plots/'];
plotname = [plotdir,sprintf('knm2_qULargeRangeScan_%s_%s.pdf',DataType,RunList)];
export_fig(f1,plotname);
fprintf('save plot to file %s \n',plotname);