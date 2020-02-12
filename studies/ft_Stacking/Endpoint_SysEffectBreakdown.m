% ------------------------------------------------------------------------------------------------%
% Switch on systematic effects one-by-one to see effect on fitted parameter
% ------------------------------------------------------------------------------------------------%
% settings:
RunList = 'StackCD100all';
SysEffect_all = {'Stat','TC','TASR','FSD','RF_EL','RF_BF','RF_RX','all'};
exclDataStart = 9;
RecomputeFlag = 'OFF';
special_label = 'Abdurashitov_ShapeOnly';%in case you want to try something, and save it
%Init Model
MRA = MultiRunAnalysis('RunList',RunList,'fixPar','1 5 6','chi2','chi2Stat',...
    'exclDataStart',exclDataStart);
belowE0 = 18575-MRA.ModelObj.qU(exclDataStart);
%Gather results
E0    = zeros(numel(SysEffect_all),1);
E0Err = zeros(numel(SysEffect_all),1);
chi2min = zeros(numel(SysEffect_all),1);
save_name = sprintf('./results/ResultE0_SystematicBreakdown_%s_%.0feV_%s.mat',RunList,belowE0,special_label);
if exist(save_name,'file')==2 && strcmp(RecomputeFlag,'OFF')
    load(save_name);
else
    %do fits
    for i=1:numel(SysEffect_all)
        SysEffect = SysEffect_all{i};
        if strcmp(SysEffect,'RF')
            SysEffects = struct('RF_EL','ON','RF_BF','ON','RF_RX','ON');
            chi2 = 'chi2CMShape';
        elseif strcmp(SysEffect,'all')
            SysEffects = struct(...
                'RF_EL','ON',...  % Response Function(RF) EnergyLoss
                'RF_BF','ON',...  % RF B-Fields
                'RF_RX','ON',...  % Column Density, inel cross ection
                'FSD','ON',...
                'TASR','ON',...
                'TCoff_RAD','ON',...
                'TCoff_OTHER','ON');
            chi2 = 'chi2CMShape';
        elseif strcmp(SysEffect,'TC')
            SysEffects = struct('TCoff_RAD','ON','TCoff_OTHER','ON');
            chi2 = 'chi2CMShape';
        elseif strcmp(SysEffect,'Stat')
            chi2 = 'chi2Stat';
        else
            SysEffects =struct(SysEffect, 'ON');
            chi2 = 'chi2CMShape';
        end
        close;
        MRA.chi2 = chi2;
        if strcmp(MRA.chi2,'chi2CMShape')
            MRA.ComputeCM('InitNormFit','ON','SysEffects',SysEffects);%,...
              %  'WGTS_TASR_RelErr',0.001);
            %        ,'ISXsection_RelErr',ISXsection_RelErr,...
            %             'WGTS_CD_MolPerCm2_RelErr',WGTS_CD_MolPerCm2_RelErr,'WGTS_B_T_RelErr',WGTS_B_T_RelErr,...
            %             'MACE_Bmax_T_RelErr',MACE_Bmax_T_RelErr,'MACE_Ba_T_RelErr',MACE_Ba_T_RelErr,...
            %             'FSDShapeES_RelErr',FSDShapeES_RelErr,'FSDShapeGS_RelErr',FSDShapeGS_RelErr,...
            %             'FSDNorm_RelErr',FSDNorm_RelErr,'WGTS_TASR_RelErr',WGTS_TASR_RelErr);
        end
        MRA.Fit('CATS','OFF');
        MRA.PlotFit('saveplot','OFF');
        E0(i)    = MRA.ModelObj.Q_i+MRA.FitResult.par(2);
        E0Err(i) = MRA.FitResult.err(2);
        chi2min(i) = MRA.FitResult.chi2min;
    end
    if exist('./results/','dir')~=7 %in case results folder doesnt exist
        mkdir ./results/
    end
    save(save_name,'E0','E0Err','chi2min','SysEffect_all');%,...
    %  'ISXsection_RelErr','WGTS_CD_MolPerCm2_RelErr','WGTS_B_T_RelErr','MACE_Bmax_T_RelErr',...
    % 'MACE_Ba_T_RelErr','FSDShapeES_RelErr','FSDShapeGS_RelErr','FSDNorm_RelErr','WGTS_TASR_RelErr');
end
%%
plotchi2 = 'OFF';

fig10 = figure(10);
switch plotchi2 
    case 'ON'
set(fig10, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 1]);
    case 'OFF'
        set(fig10, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 0.5]);
end
x = (1:numel(SysEffect_all))';
switch plotchi2
    case 'ON'
        s1 = subplot(2,1,1);
end
errorbar(x,E0,E0Err,'-o','Color',rgb('DarkCyan'),'LineWidth',3,'MarkerFaceColor',rgb('CadetBlue'),'MarkerSize',10);
ylabel('E_{0_{eff}} (eV)');
ax = gca;
ax.XTick = x;
ax.XTickLabelRotation = 15;
yticklabels(string(round(get(gca,'YTick'),5)));
xticklabels({'stat','TC','TASR','FSD','Eloss','B-fields','\rho d \sigma','combi'});
%xticklabels({'statistic','theoretical corrections','tritium activity fluctuations','final states distribution','energy loss function','magntic fields','\rhod \sigma','combi'})
xlim([0.8 numel(x)+0.2])

PrettyFigureFormat
set(gca,'FontSize',20)
grid on;
ylim([(1-1e-06)*min(E0-E0Err), (1+1e-06)*max(E0+E0Err)]);

switch plotchi2
    case 'ON'
        s2 =  subplot(2,1,2);
        plot(x,1-chi2cdf(chi2min,17),'-o','MarkerFaceColor',rgb('IndianRed'),'LineWidth',3,'MarkerSize',10,'Color',rgb('FireBrick'));
        hold on;
        plot(linspace(0,max(x)+1,numel(x)),0.05.*ones(numel(x),1),'LineStyle','--','Color',rgb('DarkSlateGray'),'LineWidth',3)
        ax = gca;
        ax.XTick = x;
        xticklabels({'stat.','TC','TASR','FSD','energy loss','B-fields','\rhod\sigma','RF','combined'})
        ylim([0 1])
        linkaxes([s1,s2],'x');
        PrettyFigureFormat;
        set(gca,'FontSize',16)
        ylabel(sprintf('p-value'));
        grid on;
end
if ~exist('./plots/png','dir')
    mkdir ./plots/png
    mkdir ./plots/pdf
end
save_name =  [sprintf('E0_SystematicBreakdown_%s_%.0feV',RunList,belowE0'),special_label];
print(fig10,['./plots/png/',save_name,'.png'],'-dpng','-r400');
%savefig(fig10,['./plots/fig/',save_name,'.fig'],'compact');
publish_figurePDF(fig10,['./plots/pdf/',save_name,'.pdf']);
