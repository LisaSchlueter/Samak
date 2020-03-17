% Run Stacking:closer look at response functions

RunList = 'KNM2_Prompt';
TwinBias_Q = 18573.70;
savedir = [getenv('SamakPath'),'knm2ana/knm2_RunStacking/results/'];
MakeDir(savedir);
savename = sprintf('%sknm2_RunStackingRF_%s_E0%.2feV.mat',savedir,RunList,TwinBias_Q);
try
    load(savename);
catch
    fprintf('File not available %s \n',savename)
    return
end

%% plot all (raw) response functions
pSingleObj = cell(nRuns,1);
qUi = 20;

f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
for i=1:nRuns
    TeMaxRaw = find((Te(i,:)-qU(i,qUi))>0.2,1);
    pSingleObj{i} = plot(Te(i,1:TeMaxRaw)-qU(i,qUi),RF(i,1:TeMaxRaw,qUi),'LineWidth',1,'Color',rgb('PowderBlue'));
    hold on;
end

%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
qUi = 20;
run = 300;
%TeMaxRaw = find((Te(run,:)-qU(run,qUi))>0.2,1);
pRaw   = plot(Te(run,:)-qU(run,qUi),RF(run,:,qUi),'LineWidth',2,'Color',rgb('PowderBlue'));
hold on;
pInter = plot(TeModel-qUModel(qUi),RFinter(run,:,qUi),'-.','LineWidth',2,'Color',rgb('Orange'));
hold off;
xlim([-0.015,0.002]);

