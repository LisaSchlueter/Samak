%% result:
NIS = 3;
ElStep2 = 0.04;
savedir = [getenv('SamakPath'),'knm2ana/knm2_ResponseFunction/results/'];
savename = sprintf('%sElossConvSanityCheck_%.0fNIS_%.2fStep.mat',savedir,NIS,ElStep2);
if exist(savename,'file')
    load(savename)
else
InitFile = @ref_FakeRun_KNM2_CD84_308x2hours;%('recomputeRF','ON');
A = RunAnalysis('RunNr',1,'FakeInitFile',InitFile,'DataType','Fake');
A.ModelObj.RFBinStep=0.01;
%A.ModelObj.InitializeRF('IntMode','Conv');

%% energy loss function
ElStep = 0.1;
E = 0:ElStep:max(A.ModelObj.Te-A.ModelObj.qUmin);
A.ModelObj.NIS = NIS;
tic;
[ElossE, ElossFunctions] = A.ModelObj.ComputeELossFunction('E',E,'ConvFlag','IntConv');
toc;
   
%%
E2 = -500:ElStep2:500;
tic;
[ElossE2, ElossFunctions2] = A.ModelObj.ComputeELossFunction('E',E2,'ConvFlag','Conv');
toc;
MakeDir(savedir);
save(savename,'ElossFunctions','ElossFunctions2','E','E2','ElStep2','ElStep');
end
%%
figure('Units','normalized','Position',[0.1,0.1,0.8,0.5])
s1 = subplot(1,2,1);
ELossFuncInt = ElossFunctions(1,:)+ElossFunctions(2,:)+ElossFunctions(3,:);
ELossFuncConv = ElossFunctions2(1,:)+ElossFunctions2(2,:).*ElStep2+ElossFunctions2(3,:).*ElStep2^2;
plot(E,ELossFuncInt,'-','LineWidth',2);
hold on;
plot(E2,ELossFuncConv,'-.','LineWidth',2);
PrettyFigureFormat('FontSize',22);
xlabel('Energy loss (eV)');
ylabel('Probability density');
leg = legend(sprintf('{\\itsimpson} integral'),sprintf('{\\itconv}() function'));
leg.EdgeColor = rgb('Silver');
hold off

s2 = subplot(1,2,2);
ELossFuncConv2 = interp1(E2,ELossFuncConv,E,'spline');
plot(E,ELossFuncInt-ELossFuncConv2,'-','LineWidth',2);
linkaxes([s1,s2],'x');
PrettyFigureFormat('FontSize',22);
xlabel('Energy loss (eV)');
ylabel('Probability density diff.');
xlim([10 45])
leg = legend(sprintf('{\\itsimpson} integral-{\\itconv}() function'));
leg.EdgeColor = rgb('Silver');

plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(gcf,plotname);
fprintf('save plot to %s \n',plotname);
