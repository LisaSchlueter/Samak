% idea: infer RW shift from data, which is not used in fit
% test uniform FPD (and ringwise later)
% fit E0 as function of upper fit boundary [-90, -40],  [-90, -30] , etc.
% goal: check if variation small and that the large range fit is realiable
% result:

RunList   = 'KNM2_RW2';%'KNM2_Prompt';
FSDFlag   = 'BlindingKNM2';
ELossFlag = 'KatrinT2';
AnaFlag   = 'StackPixel'; % uniform FPD
chi2      = 'chi2Stat';
freePar   = 'E0 Norm'; %free fit parameter
DataType  = 'Real';
RingMerge = 'Full';

RunArg = {'FSDFlag',FSDFlag,...
    'ELossFlag',ELossFlag,...
    'AnaFlag',AnaFlag,...
    'chi2',chi2,...
    'RunList',RunList,...
    'fixPar',freePar,...
    'RingMerge',RingMerge,...
    'minuitOpt','min;minos'}; 

savedir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];
savename = [savedir,sprintf('knm2_qULargeRangeScanRing_%s_%s_%s.mat',RingMerge,DataType,RunList)];

RecomputeFlag= 'OFF';
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
    
    R = RingAnalysis('RunAnaObj',T,'RingList',T.RingList);
    %%
    
    ranges       = [40,30,20 10];
    exclDataStop = zeros(numel(ranges),1);
    par          = zeros(T.nRings,T.nPar,numel(ranges));
    err          = zeros(T.nRings,T.nPar,numel(ranges));
    pval         = zeros(T.nRings,numel(ranges),1);
    BkgMean      = zeros(T.nRings,1);
    RunData      = cell(T.nRings,1);
    
    for j=1:T.nRings
        BkgMean(j) = mean(R.MultiObj(j).RunData.TBDIS(T.RunData.qU>18575)./(R.MultiObj(j).RunData.qUfrac(R.MultiObj(j).RunData.qU>18575).*R.MultiObj(j).RunData.TimeSec));
        R.MultiObj(j).ModelObj.BKG_RateSec = BkgMean(j);
        RunData{j} = R.MultiObj(j).RunData;
    end
        
    for i=1:numel(ranges)
        progressbar(i/numel(ranges));
        exclDataStop(i) =  T.GetexclDataStart(ranges(i))-1;
        
        for j=1:T.nRings
             R.MultiObj(j).exclDataStart = 1;
            R.MultiObj(j).exclDataStop =  exclDataStop(i);
        end
        
        R.FitRings;
        par(:,:,i) = R.FitResult.par;
        err(:,:,i) = R.FitResult.err;
        chi2min_tmp = R.FitResult.chi2min;
        dof_tmp     = R.FitResult.dof;
        pval(:,i)  = 1-chi2cdf(chi2min_tmp,dof_tmp);
    end
    
    upperfitrange = T.RunData.qU(exclDataStop)-T.ModelObj.Q_i;
    Q_i = T.ModelObj.Q_i;
    nRings = T.nRings;
    save(savename,'RunData','Q_i','par','err','ranges','upperfitrange','pval','exclDataStop','BkgMean','RunArg','nRings');
end
%% plot
E0    = squeeze(par(:,2,:));
E0Err = squeeze(err(:,2,:));
f1 = figure('Units','normalized','Position',[0.1,0.1,0.45,0.45]);

y1 = (par(2,:)-mean(par(2,:))).*1e3;
yErr = err(2,:).*1e3;
[l,a] = boundedline(linspace(RunData{1}.qU(11)-Q_i,0,10),zeros(10,1),800.*ones(10,1));
a.FaceColor = rgb('LightGray'); a.FaceAlpha = 0.7;
l.LineStyle = 'none';
hold on;

l1 = plot(linspace(min(upperfitrange)-2,max(upperfitrange)+2,10),zeros(10,1),'--','Color',rgb('Silver'),'LineWidth',2);
hold on
eplot = cell(nRings,1);
LineStyle = {'-o','--o',':o','-.o'};
Color = {'DodgerBlue','FireBrick','ForestGreen','GoldenRod'};
for i=1:nRings
eplot{i} = errorbar(upperfitrange,(E0(i,:)-mean(E0(i,:))).*1e3,E0Err(i,:).*1e3,LineStyle{i},'CapSize',0,...
    'LineWidth',l1.LineWidth,'Color',rgb(Color{i}),'MarkerFaceColor',rgb(Color{i}));
end
hold off;
PrettyFigureFormat('FontSize',24);
xlabel('Upper fit boundary below E_0 (eV)');
ylabel(sprintf('{\\itE}_0 - \\langle{\\itE}_0\\rangle (meV)'))
xlim([min(min(upperfitrange))-2,max(max(upperfitrange))+2])
%ylim([min(min(E0-E0Err))*1e3*1,max(max(E0+E0Err)).*1e3.*1.4])
leg = legend([eplot{:},a],'Ring 1','Ring 2','Ring 3','Ring 4',...
                    sprintf(' Data used for m_\\nu fit'));%,sprintf('\\langle{\\itE}_0\\rangle'));
leg.Location='southeast';
leg.Title.String = sprintf(' Free fit parameter: %s \n Lower fit boundary: -90 eV',freePar);
legend boxon;
leg.EdgeColor = 'none';
leg.Color = rgb('White');
leg.FontSize = 12;
% save plot
plotdir = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/plots/'];
plotname = [plotdir,sprintf('knm2_qULargeRangeScan_%s_%s.pdf',DataType,RunList)];
export_fig(f1,plotname);
fprintf('save plot to file %s \n',plotname);