RFBinStep = 0.01;
savedir = [getenv('SamakPath'),'knm2ana/knm2_ResponseFunction/results/'];
savename = sprintf('%sResponseFunctionIntegrationTest_RFBinStep%.2g.mat',savedir,RFBinStep);
if exist(savename,'file') 
    load(savename);
else
InitFile = @ref_FakeRun_KNM2_CD84_308x2hours;%('recomputeRF','ON');
A = RunAnalysis('RunNr',1,'FakeInitFile',InitFile,'DataType','Fake');
A.ModelObj.recomputeRF='ON';
A.ModelObj.RFBinStep=RFBinStep;
tic;
A.ModelObj.InitializeRF('IntMode','Conv');
toc;
RFconv = A.ModelObj.RF;
%%
RFBinStep2 = 0.1; %fixed, should matter anymore.
A.ModelObj.RFBinStep=RFBinStep2; 
tic;
A.ModelObj.InitializeRF('IntMode','Integral');
toc;
RFint = A.ModelObj.RF;
Te = A.ModelObj.Te;
qU = A.ModelObj.qU;
save(savename,'RFint','RFconv','Te','qU','InitFile','RFBinStep2')
end
%%
qUi = 27;
figure('Units','normalized','Position',[0.1,0.1,0.8,0.5])
s1 = subplot(1,2,1);
p1 = plot(Te-qU(qUi),RFconv(:,qUi),'LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
p2 = plot(Te-qU(qUi),RFint(:,qUi),'-.','LineWidth',2,'Color',rgb('Orange'));
PrettyFigureFormat('FontSize',24);
xlabel('Surplus energy (eV)');
ylabel('Transmission probability');
leg = legend([p2,p1],sprintf('{\\itsimpson} integral'),sprintf('{\\itconv}() function'),...
    'Location','northwest');
leg.EdgeColor = rgb('Silver');
leg.Title.String = sprintf('{\\itqU}=%.1f eV',qU(qUi));
leg.Title.FontWeight = 'normal';
title(sprintf('Bin width for {\\itconv}() = %.3f',RFBinStep),'FontWeight','normal');

s2 = subplot(1,2,2);
plot(Te-qU(qUi),abs(RFconv(:,qUi)-RFint(:,qUi)),'LineWidth',2);
hold on;
plot(Te-mean(qU),abs(RFconv(:,:)-RFint(:,:)),':','LineWidth',1);
PrettyFigureFormat('FontSize',24);
xlabel('Surplus energy (eV)');
ylabel('Abs. transmission probability diff.')
leg = legend(sprintf('{\\itsimpson}-{\\itconv}'),'Location','northwest');
leg.EdgeColor = rgb('Silver');
linkaxes([s1,s2],'x');
xlim([-2 40])
ylim([0 1e-3])
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(gcf,plotname);
fprintf('save plot to %s \n',plotname);
