

savedir = sprintf('%sknm12Combi/knm12_General/results/',getenv('SamakPath'));
savename = sprintf('%sknm12_DataModel.mat',savedir);

if exist(savename,'file')
    load(savename);
else
    fprintf('run knm12_Signal2Background \n')
    return
end

if ~exist('TBDIS_i1','var')
    
    %% init models
    Model1.BKG_RateSec_i = Model1.BKG_RateSec;
    Model1.normFit_i     = Model1.normFit;
    Model1.Q_i           = Model1.Q;
    qU1                  = Model1.qU;
    TBDIS_i1             = Model1.TBDIS;
    
    Model2.BKG_RateSec_i = Model2.BKG_RateSec;
    Model2.normFit_i     = Model2.normFit;
    Model2.Q_i           = Model2.Q;
    qU2                  = Model2.qU;
    TBDIS_i2             = Model2.TBDIS;
    
    %% calculate neutrino mass signal for some nu-masses
    mNuSq_i = [0 0.1 0.5 1];
    
    TBDIS1 = zeros(numel(qU1),numel(mNuSq_i));
    TBDIS2 = zeros(numel(qU2),numel(mNuSq_i));
    
    for i=1:numel(mNuSq_i)
        Model1.mnuSq_i = mNuSq_i(i);
        Model1.ComputeTBDDS;
        Model1.ComputeTBDIS;
        TBDIS1(:,i) =  Model1.TBDIS;
        
        Model2.mnuSq_i = mNuSq_i(i);
        Model2.ComputeTBDDS;
        Model2.ComputeTBDIS;
        TBDIS2(:,i) =  Model2.TBDIS;
    end
    save(savename,'TBDIS_i1','TBDIS_i2','TBDIS1','TBDIS2','qU1','qU2','mNuSq_i','-append')
end

%% plot
Colors = [rgb('ForestGreen');rgb('DodgerBlue');rgb('Orange');rgb('FireBrick')];
FontSize = 24;

f1 = figure('Units','normalized','Position',[1.1,1.1,0.5,0.7]);

%% signal-to-background 
% select range [E0-40,40] eV
StartIdx1 = 13;
StartIdx2 = 11;
StopIdx1  = 34;
StopIdx2  = 33;
% absolute (sum) model
SigSum2 = sum(SignalCount2(StartIdx2:StopIdx2));
BkgSum2 = sum(BkgCounts2(StartIdx2:StopIdx2));
Sig2Bkgtot2 = SigSum2/BkgSum2;

SigSum1 = sum(SignalCount1(StartIdx1:StopIdx1));
BkgSum1 = sum(BkgCounts1(StartIdx1:StopIdx1));
Sig2Bkgtot1 = SigSum1/BkgSum1;

% plot
s1 = subplot(4,1,[1:2]);
plot(qU1-18574,ones(39,1),':','LineWidth',2.5,'Color',rgb('Silver'));
hold on
qUinter = linspace(18573-7-50,18573.7-2.8,1e2);
Sig2Bkg1_inter = interp1(qU1,Sig2Bkg1,qUinter,'spline');
Sig2Bkg2_inter = interp1(qU2,Sig2Bkg2,qUinter,'spline');
p1 = plot(qUinter-Model1.Q,Sig2Bkg1_inter,'-','LineWidth',3,'Color',rgb('ForestGreen'),'MarkerSize',18);
p2 = plot(qUinter-Model2.Q,Sig2Bkg2_inter,'-.','LineWidth',3,'Color',p1.Color,'MarkerSize',18);


ylabel('Signal-to-background');
xlim([-40 0]);
xticklabels('');
ylim([1e-03 1e3]);
yticks([1e-2 1 1e2])
PrettyFigureFormat('FontSize',FontSize);

set(gca,'YScale','log');
ax1 = gca;

leg = legend([p1,p2],'KNM1','KNM2','Location','southwest');
legend boxoff
leg.Title.String = sprintf('{\\itm}_\\nu^2 = 0 eV^2');
leg.Title.FontWeight = 'normal'; leg.Title.FontSize = leg.FontSize;
t1 = text(-39.5,550,'a)','FontName',get(gca,'FontName'),'FontSize',get(gca,'FontSize'));

% plot 2
s2 = subplot(4,1,[3:4]);

legH = cell(numel(mNuSq_i)-1,1);
legStr = cell(numel(mNuSq_i)-1,1);

for i=1:numel(mNuSq_i)-1
    if i==1
        qU_inter1 = linspace(min(qU1)-Model1.Q,0,1e3);  
        TBDIS1_tmp =  interp1(qU1-Model1.Q,TBDIS1(:,i)./TBDIS_i1,qU_inter1,'lin');
        
         qU_inter2 = linspace(min(qU2)-Model2.Q,0,1e3);  
        TBDIS2_tmp =  interp1(qU2-Model2.Q,TBDIS2(:,i)./TBDIS_i2,qU_inter2,'lin');
    else
        qU_inter1 = linspace(min(qU1)-Model1.Q,-2.38,1e3);
        TBDIS1_tmp =  interp1(qU1-Model1.Q,TBDIS1(:,i)./TBDIS_i1,qU_inter1,'spline');
        
         qU_inter2 = linspace(min(qU2)-Model2.Q,0,1e3);
        TBDIS2_tmp =  interp1(qU2-Model2.Q,TBDIS2(:,i)./TBDIS_i2,qU_inter2,'spline');
    end
   
    legH{i} = plot(qU_inter1,TBDIS1_tmp,'LineWidth',2,'Color',Colors(i,:));
    legStr{i} = sprintf('{\\itm}_\\nu^2 = %.2g eV^2',mNuSq_i(i));
    hold on;  
    plot(qU_inter2,TBDIS2_tmp,'-.','LineWidth',2,'Color',Colors(i,:));
end

%%

xlim([-40 0]);

ylim([0.989 1.0015])
ylabel('Neutrino-mass signal');
xlabel(sprintf('Retarding energy - {\\itE}_0 (eV)'));
ax2 = gca;
legend([legH{:}],legStr,'Location','southwest');
legend box off
PrettyFigureFormat('FontSize',FontSize);
t2 = text(t1.Position(1),1.00082,'b)','FontName',get(gca,'FontName'),'FontSize',get(gca,'FontSize'));

linkaxes([s1,s2],'x');
% position
ax1.Position(2) = 0.54;
ax1.Position(4) = 0.42;
ax2.Position(4) =0.42;
ax1.Position(1) = 0.2;
ax2.Position(1) = ax1.Position(1);
ax1.YLabel.Position(1) = ax2.YLabel.Position(1);
%%
fprintf('KNM1 signal to background ratio of 2 at %.1f eV \n',interp1(Sig2Bkg1(25:31),qU1(25:31)-Model1.Q,1,'spline'))
fprintf('KNM2 signal to background ratio of 2 at %.1f eV \n',interp1(Sig2Bkg2(24:28),qU2(24:28)-Model2.Q,1,'spline'))

%%
pltdir = [getenv('SamakPath'),'knm12Combi/knm12_General/plots/'];
MakeDir(pltdir);
pltfile = sprintf('%sknm12_NuMassSignal.pdf',pltdir);
export_fig(pltfile);
fprintf('Save plot to %s \n',pltfile)