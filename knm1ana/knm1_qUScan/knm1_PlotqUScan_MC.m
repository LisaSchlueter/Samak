% plot knm1 qU scan with MC expectation

qURange  = [95,20];
chi2 = 'chi2CMShape';
nSamples = 1000;

savedir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/results/'];

% load regular qU Scan
savename = sprintf('%sknm1_qUScan_Mini_%.0feV_to_%.0feV_%s_NP1.064.mat',...
    savedir,qURange(1),qURange(2),chi2);
d = importdata(savename);

% load randomized qU Scan
savenameMC = sprintf('%sknm1_qUScan_MC_%.0feV_to_%.0feV_%s_NP%2g_%.0fsamples.mat',...
    savedir,qURange(1),qURange(2),chi2,1.064,nSamples);
dMC = importdata(savenameMC);

%% display
close all;

fig12345 = figure('Renderer','painters');
set(fig12345, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.3]);

[l,a] = boundedline(dMC.qU(dMC.exclDataStart_v)-18574,dMC.mNuSq_mu+d.parqU(1,13),dMC.mNuSq_std);
hold on;
l.delete;
a.FaceColor = rgb('LightGray');
a.FaceAlpha = 0.8;

[parqU, errqU, chi2qU, dofqU,e1,pref] =  d.M.qUScan('qURange',qURange,...
    'saveplot','OFF',...
    'RecomputeFlag','OFF',...
    'CorrMean','OFF',...
    'HoldOn','ON',...
    'RelFlag','OFF',...
    'RefLine',40,...
    'PlotResult','ON',...
    'saveStr','_KatrinT2');
e1.Color = rgb('DodgerBlue');

%%
leg = legend([pref,a],'Standard analysis range',...
    sprintf('Expected 1\\sigma variation'));
leg.ItemTokenSize = [30,13];
leg.Location = 'southwest';


%% correlation matrix
%close all
fig12345 = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.4, 0.5]);


FontSize = 22;
PltPar = 1;


switch PltPar
    case 1
        mNuSq_abs = dMC.mNuSq_abs';
        mNuSq_abs(dMC.Idx_rm) = NaN;
        CorrMat = corr(mNuSq_abs,'rows','complete');
        tStr = sprintf('{\\itm}_\\nu^2 ');
    case 2
        E0_abs = dMC.E0_abs';
        E0_abs(dMC.Idx_rm) = NaN;
        CorrMat = corr(E0_abs,'rows','complete');
        tStr = sprintf('{\\itE}_0^{fit} ');
    case 3
        B_abs = dMC.B_abs';
        B_abs(dMC.Idx_rm) = NaN;
        CorrMat = corr(B_abs,'rows','complete');
        tStr = sprintf('{\\itB}_{base} ');
    case 4
        N_abs = dMC.N_abs';
        N_abs(dMC.Idx_rm) = NaN;
        CorrMat = corr(N_abs,'rows','complete');
        tStr = sprintf('{\\itN}_{sig.}');
end

imagesc(CorrMat);
pbaspect([1 1 1])

% xy ticks
TickLabels = strings(23,1);
for i=1:23
    if i==23 || i==1 || i==13%mod(i,6)==1 ||
        TickLabels{i} = sprintf('%.0f eV',dMC.qU(dMC.exclDataStart_v(i))-18574);
    end
end
xticks([1:23]);
yticks([1:23]);
xticklabels(TickLabels);
yticklabels(TickLabels);
xlabel('Lower fit boundary (bin number)')
ylabel('Lower fit boundary (bin number)')

PrettyFigureFormat('FontSize',FontSize);
set(gca,'XMinorTick','off');
set(gca,'YMinorTick','off');
ax = gca;

cb = colorbar;
cb.Label.String = sprintf('Correlation coefficient (%s)',tStr);
cb.Label.FontSize = ax.XLabel.FontSize;
%
pltdir = [getenv('SamakPath'),'knm1ana/knm1_qUScan/plots/'];
MakeDir(pltdir);
pltname = sprintf('%sknm1_qUScan_CorrMat_Par%.0f.pdf',pltdir,PltPar);
export_fig(pltname);
fprintf('save plot %s \n',pltname);



