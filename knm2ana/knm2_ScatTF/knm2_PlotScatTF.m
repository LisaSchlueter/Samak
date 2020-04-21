% plot transmission function with different options

%% settings
RecomputeFlag = 'OFF';
%% load if possible
savedir = [getenv('SamakPath'),'knm2ana/knm2_ScatTF/results/'];
%% OFF - OFF
AngularTFFlag = 'OFF';
SynchrotronFlag = 'OFF';
savename1 = sprintf('%sknm2_TF_Sync%s_Scat%s.mat',savedir,SynchrotronFlag,AngularTFFlag);
if exist(savename1,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename1)
else
    A = ref_FakeRun_KNM2_RFcomparison('SynchrotronFlag',SynchrotronFlag,...
                                       'AngularTFFlag',AngularTFFlag);
    qu = 18545;
    te = (0:0.01:3)+qu;
    tic;
    MaceTF_SyncOFF_ScatOFF = A.ComputeMaceTF(te',qu);
    timeOFFOFF  = toc;
    save(savename1,'te','qu','MaceTF_SyncOFF_ScatOFF','timeOFFOFF');
end

%% ON - OFF
AngularTFFlag = 'ON';
SynchrotronFlag = 'OFF';
savename3 = sprintf('%sknm2_TF_Sync%s_Scat%s.mat',savedir,SynchrotronFlag,AngularTFFlag);
if exist(savename3,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename3)
else
    A = ref_FakeRun_KNM2_RFcomparison('SynchrotronFlag',SynchrotronFlag,...
                                       'AngularTFFlag',AngularTFFlag);
    qu = 18545;
    te = (0:0.01:3)+qu;
    tic;
    MaceTF_SyncOFF_ScatON = A.ComputeMaceTF(te',qu);
    timeScatON = toc;
    save(savename3,'te','qu','MaceTF_SyncOFF_ScatON','timeScatON');
end

%% OFF - ON
AngularTFFlag = 'OFF';
SynchrotronFlag = 'ON';
savename3 = sprintf('%sknm2_TF_Sync%s_Scat%s.mat',savedir,SynchrotronFlag,AngularTFFlag);
if exist(savename3,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename3)
else
    A = ref_FakeRun_KNM2_RFcomparison('SynchrotronFlag',SynchrotronFlag,...
                                       'AngularTFFlag',AngularTFFlag);
    qu = 18545;
    te = (0:0.01:3)+qu;
    tic;
    MaceTF_SyncON_ScatOFF= A.ComputeMaceTF(te',qu);
    timeSyncON = toc;
    save(savename3,'te','qu','MaceTF_SyncON_ScatOFF','timeSyncON');
end
%% ON - ON
AngularTFFlag = 'ON';
SynchrotronFlag = 'ON';
savename4 = sprintf('%sknm2_TF_Sync%s_Scat%s.mat',savedir,SynchrotronFlag,AngularTFFlag);
if exist(savename4,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename4)
else
    A = ref_FakeRun_KNM2_RFcomparison('SynchrotronFlag',SynchrotronFlag,...
                                       'AngularTFFlag',AngularTFFlag);
    qu = 18545;
    te = (0:0.01:3)+qu;
    tic;
    MaceTF_SyncON_ScatON = A.ComputeMaceTF(te',qu);
    timeSyncONScatON = toc;
    save(savename4,'te','qu','MaceTF_SyncON_ScatON','timeSyncONScatON');
end

%% plot anuglar scattering ON - OFF
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
LineWidth = 2;
nScat = 1;
Sync = 'ON';
if strcmp(Sync,'ON')
    p1 = plot(te-qu,MaceTF_SyncON_ScatOFF,'-','LineWidth',LineWidth,'Color',rgb('Orange'));
    hold on;
    p2 = plot(te-qu,MaceTF_SyncON_ScatON(nScat,:),'-.','LineWidth',LineWidth,'Color',rgb('DodgerBlue'));
else
    p1 = plot(te-qu,MaceTF_SyncOFF_ScatOFF,'-','LineWidth',LineWidth,'Color',rgb('Orange'));
    hold on;
    p2 = plot(te-qu,MaceTF_SyncOFF_ScatON(nScat,:),'-.','LineWidth',LineWidth,'Color',rgb('DodgerBlue'));
end
leg = legend('Isotropic','Angular scattering');
leg.EdgeColor = rgb('Silver');
leg.Title.String = sprintf('Synchrotron  %s',Sync);
leg.Title.FontWeight = 'normal';
leg.Location = 'northwest';
PrettyFigureFormat('FontSize',24);
xlabel(sprintf('{\\itE - qU} (eV)'));
ylabel('Transmission probability');
ylim([-0.05 1.05])
plotname = sprintf('%sScatTransmission_Synchr%s.png',savedir,Sync);
print(gcf,plotname,'-dpng','-r400');
fprintf('save plot to %s \n',plotname);

