% look at variations of column density within a scan
% plot for PhD thesis
savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    return
end

qU = R.RunData.qU(R.exclDataStart:end,:);
CD = R.SingleRunData.WGTS_CD_MolPerCm2_SubRun(R.exclDataStart:end,:);
CD_scatter = reshape((CD-mean(CD,2))./mean(CD,2),1,R.nRuns.*numel(qU))';
qU_scatter = repmat(qU,R.nRuns,1);
%%
GetFigure

sp = dscatter(qU_scatter,1e2*(CD_scatter),'Msize',20); % density scatter plot
hold on;
%plot(qU,1e2.*(mean(CD,2)-R.ModelObj.WGTS_CD_MolPerCm2)./R.ModelObj.WGTS_CD_MolPerCm2,'.r','MarkerSize',12);
plot(qU,1e2*(+std(CD,0,2))./mean(CD,2),':r','LineWidth',2);
plot(qU,1e2*(-std(CD,0,2))./mean(CD,2),':r','LineWidth',2);
plot(qU,1e2*zeros(size(mean(CD,2))),'-r','LineWidth',2);
xlabel('Retarding energy');
ylabel('Relative column density');%sprintf('[\\rhod - \\langle\\rhod\\rangle]/\\langle\\rhod\\rangle (%%)'));
PrettyFigureFormat('FontSize',22);
c = colorbar;
c.Label.String = 'Scan density'; c.Label.FontSize = ax.XLabel.FontSize;
colormap(flipud(colormap('winter')));
ax = gca;
ax.XAxis.Exponent = 0;
xlim([min(qU)-3 max(qU)+3]);
xticks(18540:20:18620);


%%
MeanCD_StackedRun = R.ModelObj.WGTS_CD_MolPerCm2;%   .mean(mean(CD,2));
MeanCD_SingleRun  = mean(CD);
CDcorr = CD.*MeanCD_StackedRun./MeanCD_SingleRun; 

WGTS_TASR_AbsErr = std(CDcorr,0,2)./sqrt(R.nRuns);%std(CDcorr,0,2)./sqrt(numel(obj.StackedRuns));
AbsErr = sqrt(cov(CDcorr'))./sqrt(R.nRuns); 
RelErr = WGTS_TASR_AbsErr./MeanCD_StackedRun;
                        % relative uncertainty -> fractional covariance
                        % WGTS_TASR_RelErr = WGTS_TASR_AbsErr./mean(SubRunActivity,2);

