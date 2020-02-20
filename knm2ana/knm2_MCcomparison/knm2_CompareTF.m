% compare tranmission functions
SynchrotronFlag = 'OFF';

% Samak TF
savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
switch SynchrotronFlag
    case 'ON'
        savename = [savedir,'Samak_TF_Synchr.mat'];
    case 'OFF'
        savename = [savedir,'Samak_TF_NoSynchr.mat'];
end

if exist(savename,'file') 
    load(savename);
else
    A = ref_FakeRun_KNM2_RFcomparison('SynchrotronFlag',SynchrotronFlag);
    qu = 18545;
    te = (0:0.01:3)+qu;
    MaceTF = A.ComputeMaceTF(te,qu);
    %plot(te-qu,MaceTF);
    Write2Txt('filename',strrep(savename,'.mat',''),'nCol',2,'variable',[te;MaceTF]);
    save(savename,'MaceTF','qu','te','A');
end
xS = (te-qu)';

% load fitrium
dF = importdata([savedir,'TF_Fitrium_synchrotron.dat']);
xF = dF.data(:,1); % surplus energy
switch SynchrotronFlag
    case 'ON'
        yF = dF.data(:,2); % transmission with synchrotron
    case 'OFF'
        yF = dF.data(:,3); % transmission without synchrotron
end

% load Kafit
dK = importdata([savedir,'TF_Kafit_synchrotron.dat']);
xK = dK(:,1); % surplus energy
yK = dK(:,2); % transmission with synchrotron
Index = ismember(round(xK,5),round(xS,5));
xK = xK(Index);
yK = yK(Index);

%% plot difference
f1 = figure('Units','normalized','Position',[0.5,0.1,0.5,0.5]);
l = plot(linspace(-5,90,100),zeros(100,1),'-','LineWidth',2,'Color',rgb('Black'));
hold on;
pF = plot(xS,MaceTF'-yF,'-.','LineWidth',2.5,'Color',rgb('DodgerBlue'));
if strcmp(SynchrotronFlag,'ON')
    pK = plot(xS,MaceTF'-yK,'-','LineWidth',2.5,'Color',rgb('GoldenRod'));
    pFK = plot(xS,yF-yK,':','LineWidth',2.5,'Color',rgb('IndianRed'));
end
xlabel(sprintf('Energy - %.0f (eV)',qu));
ylabel('Probability diff.');
PrettyFigureFormat('FontSize',22);
if strcmp(SynchrotronFlag,'ON')
    leg = legend([pK,pF,pFK],'Samak - KaFit','Samak - Fitrium','Fitrium - KaFit');
    leg.Title.String = 'Synchrotron ON';
    leg.Location = 'northwest';
else
    leg = legend(pF,'Samak - Fitrium');
    leg.Title.String = 'Synchrotron OFF';
    leg.Location = 'southwest';
end
leg.EdgeColor = rgb('Silver');
xlim([min(xS),max(xS)]);

plotdir = strrep(savedir,'results','plots');
savename = sprintf('%sTF_Diff_Synchrotron%s',plotdir,SynchrotronFlag);
export_fig(f1,[savename,'.pdf']);
print(f1,[savename,'.png'],'-dpng','-r300');



