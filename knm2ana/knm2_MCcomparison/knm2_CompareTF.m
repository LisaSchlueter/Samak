% compare tranmission functions
SynchrotronFlag = 'ON';
RecomputeFlag = 'OFF';
AngTF = 'ON';
SaveTxt = 'ON';

% Samak TF
savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/Transmission/'];
savename = sprintf('%sSamakTF_Sync%s_ScatTF%s.mat',savedir,SynchrotronFlag,AngTF);
 
if strcmp(SynchrotronFlag,'OFF') && strcmp(AngTF,'ON')
    fprintf(2,'angular tf only with syncrotron from fitrium \n')
    return
end
if exist(savename,'file')  && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    A = ref_FakeRun_KNM2_RFcomparison('SynchrotronFlag',SynchrotronFlag,...
                                       'AngularTFFlag',AngTF);
    qu = 18545;
    te = (0:0.01:3)+qu;
    MaceTF = A.ComputeMaceTF(te',qu);
    
    save(savename,'MaceTF','qu','te','A');
    if strcmp(SaveTxt,'ON')
        if strcmp(AngTF,'ON')
            Write2Txt('filename',strrep(savename,'.mat',''),'nCol',3,'variable',[te;MaceTF(1:2,:)]);
        else
            Write2Txt('filename',strrep(savename,'.mat',''),'nCol',2,'variable',[te;MaceTF']);
        end
    end
end
xS = (te-qu)';

if strcmp(AngTF,'ON')
    MaceTF1 = MaceTF(2,:)'; % one scattering
    MaceTF2 = MaceTF(3,:)'; % two scattering
    MaceTF  = MaceTF(1,:)'; % zero scattering
end

% load fitrium
if strcmp(AngTF,'ON') || strcmp(SynchrotronFlag,'ON')
    dF = importdata([savedir,'TF_Fitrium_sync_scat.dat']); % energy, no angular, 0 scat, 1 scat,
    xF = dF(:,1); % surplus energy
    Index = ismember(round(xF,5),round(xS,5));
    xF = xF(Index);
    
    if strcmp(AngTF,'ON')
        yF = dF(:,3); % 0 scat
        yF1 = dF(:,4); % 1 scat
        yF1 = yF1(Index);
    else
        yF = dF(:,2);
    end
    yF = yF(Index);
else
    dF = importdata([savedir,'TF_Fitrium_sync.dat']);
    switch SynchrotronFlag
        case 'ON'
            yF = dF.data(:,2); % transmission with synchrotron
        case 'OFF'
            yF = dF.data(:,3); % transmission without synchrotron
    end
    xF = dF.data(:,1); % surplus energy
end


%% load Kafit
if strcmp(SynchrotronFlag,'ON') && strcmp(AngTF,'OFF')
    dK = importdata([savedir,'TF_Kafit_sync.dat']);
elseif strcmp(SynchrotronFlag,'OFF') && strcmp(AngTF,'OFF')
    dK = importdata([savedir,'TF_Kafit.dat']);
elseif strcmp(SynchrotronFlag,'ON') && strcmp(AngTF,'ON')
    dK = importdata([savedir,'TF_Kafit_sync_scat.dat']); 
end
xK = dK(:,1); % surplus energy
yK = dK(:,2); % transmission with synchrotron
Index = ismember(round(xK,5),round(xS,5));
xK = xK(Index);
yK = yK(Index);

%% plot difference
f1 = figure('Units','normalized','Position',[0.5,0.1,0.5,0.48]);
l = plot(linspace(-5,90,100),zeros(100,1),'-','LineWidth',2,'Color',rgb('Black'));
hold on;
pF = plot(xS,MaceTF-yF,'-','LineWidth',2.5,'Color',rgb('DodgerBlue'));

if  strcmp(AngTF,'OFF')
    pK = plot(xS,MaceTF-yK,'--','LineWidth',2.5,'Color',rgb('GoldenRod'));
    pFK = plot(xS,yF-yK,':','LineWidth',2.5,'Color',rgb('IndianRed'));
else
     hold on;
     pF1 = plot(xF,MaceTF1-yF1,'-.','Color',rgb('Orange'),'LineWidth',2.5);
end

xlabel(sprintf('Surplus energy (eV)'));
ylabel('Probability diff.');
PrettyFigureFormat('FontSize',22);
if strcmp(AngTF,'OFF')
    leg = legend([pK,pF,pFK],'Samak - KaFit','Samak - Fitrium','Fitrium - KaFit');
    AngStr = 'isotropic transmission';
else
    leg = legend([pF,pF1],'Samak - Fitrium (0 scat.)','Samak - Fitrium (1 scat.)');
    AngStr = 'non-isotropic transmission';
end
if strcmp(SynchrotronFlag,'ON')
    SyncStr = 'Synchrotron radiation';
else
    SyncStr = 'Without synchrotron radiation';
end

 title(sprintf('%s and %s',SyncStr,AngStr),...
       'FontWeight','normal','FontSize',get(gca,'FontSize'));
leg.Title.FontWeight = 'normal';
leg.EdgeColor = rgb('Silver');
leg.Location = 'northwest';
xlim([min(xS),max(xS)]);
ylim([-5 1.5].*1e-07);
set(gca,'XMinorTick','off');

plotdir = strrep(savedir,'results/Transmission','plots');
MakeDir(plotdir);
savename = sprintf('%sTF_Diff_Sync%s_ScatTF%s',plotdir,SynchrotronFlag,AngTF);
export_fig(f1,[savename,'.pdf']);
print(f1,[savename,'.png'],'-dpng','-r300');
fprintf('save plot to %s \n',savename)



