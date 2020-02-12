RunList = 'KNM1';
DataType = 'Real';
chi2 = 'chi2CMShape';
SysBudet = 22;
fixPar = '5 6 7 8 9 10 11';%5 6 7 8 9 10 11';%'11'
firstPoint = 2;
lastPoint = 20;
FSDFlag ='SibilleFull';
M = MultiRunAnalysis('RunList',RunList,...
    'chi2','chi2Stat','DataType',DataType,...
    'fixPar',fixPar,...
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1.064,...
    'minuitOpt','min ; migrad',...
    'FSDFlag',FSDFlag,...
    'ELossFlag','KatrinT2',...
    'SysBudget',22);
M.chi2 = chi2;
%%
if strcmp(fixPar,'5 6 7 8 11')
    M.ComputeCM('FSDNorm_RelErr',0);
else
    M.ComputeCM;
end
%%
[parqU, errqU, chi2qU, dofqU] = ...
    M.qUScan('firstPoint',firstPoint,'lastPoint',lastPoint,'saveplot','OFF',...
    'RecomputeFlag','OFF','CorrMean','OFF','HoldOn','OFF','RelFlag','OFF');
% sibille
%%
PrettyFigureFormat;
set(gca,'FontSize',24)
leg = legend('Sibille','Sibille + FSD Normalization free');%'Saenz')
leg.Location='southwest';
leg.FontSize = 24;

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4)- ti(2) - ti(4)-0.01;
ax.Position = [left+0.15 bottom ax_width-0.4 ax_height-0.1];

%%
savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/plots/'];
savename = [savedir,sprintf('qUScan_%s_%s_%s_%s_SibilleFSDfree.pdf',...
M.RunData.RunName,M.DataType,M.chi2,'567891011')];

publish_figurePDF(gcf,savename);

%% additional plot: TT ground state proability
fig12345 = figure('Renderer','painters');
set(fig12345, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.6]);
ColorArg = {'MarkerFaceColor',rgb('SlateGray'),'Color',rgb('DarkSlateGray')};

Mode = 'T2';
switch Mode
    case 'T2'
        y = flip(parqU(9,:))+M.ModelObj.TTNormGS_i;%
        yErr = flip(errqU(9,:));
        yref = M.ModelObj.TTNormGS_i;
        str = sprintf('T_2');
    case 'DT'
        y = flip(parqU(5,:))+M.ModelObj.TTNormGS_i;%
        yErr = flip(errqU(5,:));
        yref = M.ModelObj.DTNormGS_i;
        str = 'DT';
    case 'HT'
        y = flip(parqU(7,:))+M.ModelObj.TTNormGS_i;%
        yErr = flip(errqU(7,:));
        yref = M.ModelObj.HTNormGS_i;
        str = 'HT';
end

x =flip(M.RunData.qU(2:20,1))-18575;%obj.ModelObj.Q_i;

[l1,a1] = boundedline(linspace(min(x)-5,max(x)+5,numel(x)),yref.*ones(numel(x),1),0.01.*M.ModelObj.TTNormGS_i.*ones(numel(x),1));
l1.LineWidth=2;
hold on;
e1 = errorbar(x,y,yErr,'o',...
    ColorArg{:},'LineWidth',2.5,'MarkerSize',12);
e1.CapSize = 0;
xlabel(sprintf('Lower fit boundary below E_0 (eV)'));
PrettyFigureFormat;
ylabel(sprintf('Ground state probability %s',str))

leg = legend([e1,l1,a1],'Fit result','Expectation (Sibille FSD)',sprintf('1 %% uncertainty band'));
legend boxoff; leg.Location='northwest';
%ylim([0.563,0.59]);

if strcmp(M.DataSet,'Knm1')
    savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/plots/'];
elseif contains(obj.DataSet,'FirstTritium')
    savedir = [getenv('SamakPath'),'first-tritium-studies/ft_qUScan/plots/'];
    MakeDir(savedir);
else
    savedir = './plots/';
end
     
savename = [savedir,sprintf('%s_GSProb_qUScan.pdf',Mode)];
publish_figurePDF(fig12345,savename);

                        
                        


