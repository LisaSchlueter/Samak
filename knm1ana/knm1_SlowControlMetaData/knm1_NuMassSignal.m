


savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_NuMassSignal.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    range = 40; % eV below the endpoint
    %% create MultiRunAnalysis object
    RunAnaArg =  {'RunList','KNM1',... % runlist defines which runs are analysed -> set MultiRunAnalysis.m -> function: GetRunList()
        'chi2','chi2Stat',...                 % uncertainties: statistical or stat + systematic uncertainties
        'DataType','Twin',...                 % can be 'Real' or 'Twin' -> Monte Carlo
        'fixPar','E0 Norm Bkg',...        % free Parameter!!
        'NonPoissonScaleFactor',1,...        % background uncertainty are (not) enhanced
        'minuitOpt','min ; minos',...         % technical fitting options (minuit)
        'FSDFlag','SibilleFull',...          % final state distribution
        'ELossFlag','KatrinT2',...            % energy loss function
        'SysBudget',22,...                    % defines syst. uncertainties -> in GetSysErr.m;
        'RadiativeFlag','ON',...              % theoretical radiative corrections applied in model
        'DopplerEffectFlag','OFF',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'DopplerEffectFlag','OFF'};   % already included inside FSD for knm1!!!!
    %%
    T = MultiRunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range); % set region of interest
    
    %% re-set init values to match Twin data
    T.Fit;
    T.ModelObj.BKG_RateSec_i = T.ModelObj.BKG_RateSec;
    T.ModelObj.normFit_i = T.ModelObj.normFit;
    T.ModelObj.Q_i = T.ModelObj.Q; 
    
    T.exclDataStart = T.exclDataStart-1;
    TBDIS_i = T.ModelObj.TBDIS(T.exclDataStart:end);  
     %%
   
    mNuSq_i = [0 0.1 0.5 1];
    qU    = T.ModelObj.qU(T.exclDataStart:end);
    Te    = T.ModelObj.Te;
    %%
    TBDIS = zeros(numel(qU),numel(mNuSq_i));
    TBDDS = zeros(numel(Te),numel(mNuSq_i));
    
    for i=1:numel(mNuSq_i)
        T.ModelObj.mnuSq_i = mNuSq_i(i);
        T.ModelObj.ComputeTBDDS;
        TBDDS(:,i) = T.ModelObj.TBDDS;
        T.ModelObj.ComputeTBDIS;
        TBDIS(:,i) = T.ModelObj.TBDIS(T.exclDataStart:end);
    end
  %%  
    TBDIS_Bkg = T.ModelObj.BKG_RateSec.*T.ModelObj.qUfrac(T.exclDataStart:end).*T.ModelObj.TimeSec;
    TBDIS_Sig = TBDIS- TBDIS_Bkg;
    TBDIS_i_Sig = TBDIS_i- TBDIS_Bkg;
    
    save(savefile,'T','TBDIS','TBDIS_i',...
        'TBDIS_Bkg','TBDIS_Sig','TBDIS_i_Sig',...
        'mNuSq_i','qU','Te','TBDDS');
end

%%
close all;
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.7]);
s1 = subplot(4,1,[1:3]);

Colors = [rgb('ForestGreen');rgb('DodgerBlue');rgb('Orange');rgb('FireBrick')];

legH = cell(numel(mNuSq_i),1);
legStr = cell(numel(mNuSq_i),1);

for i=1:numel(mNuSq_i)
    if i==1
        qU_inter = linspace(min(qU)-T.ModelObj.Q,0,1e3);  
         TBDIS_tmp =  interp1(qU-T.ModelObj.Q,TBDIS(:,i)./TBDIS_i,qU_inter,'lin');
    else
        qU_inter = linspace(min(qU)-T.ModelObj.Q,-2.38,1e3);
         TBDIS_tmp =  interp1(qU-T.ModelObj.Q,TBDIS(:,i)./TBDIS_i,qU_inter,'spline');
    end
   
    legH{i} = plot(qU_inter,TBDIS_tmp,'LineWidth',2,'Color',Colors(i,:));
    legStr{i} = sprintf('{\\itm}_\\nu^2 = %.2g eV^2',mNuSq_i(i));
    hold on;  
end

PrettyFigureFormat('FontSize',24);
xlim([-40 0]);
ylim([0.989,1.001])
ylabel('Neutrino-mass signal');
xticklabels('');
ax1 = gca;
leg = legend([legH{:}],legStr,'Location','southwest');
legend box off
%%
t = text(-39.5,1.0005,'a)','FontName',get(gca,'FontName'),'FontSize',get(gca,'FontSize'));
%%
s2 = subplot(4,1,4);
b1 = bar(qU-T.ModelObj.Q,TBDIS_i_Sig./TBDIS_Bkg,'FaceColor',rgb('Silver'),'EdgeColor',rgb('Silver'));
hold on;
hold on;
plot(linspace(-40,0,10),ones(10,1),':k','LineWidth',2)
pNone = plot(NaN,NaN,'.w');
xlabel('Retarding energy - 18573.7 (eV)');
ylabel(sprintf('{\\itR}_{sig.}/{\\itR}_{bkg.}'));
set(gca,'YScale','log');
PrettyFigureFormat('FontSize',24);
ax2 = gca;

xlim([-40 0]);
ylim([1e-04 700]);
linkaxes([s1,s2],'x');
leg = legend(pNone,sprintf('Signal-to-background ratio for {\\itm}_\\nu^2 = 0 eV^2'),...
    'Location','northeast');
leg.ItemTokenSize(2) = 10;
legend box off
% position
ax1.Position(2) = 0.415;
ax1.Position(4) = 0.58;
ax2.Position(4) =ax1.Position(4)/2;
t = text(-39.5,220,'b)','FontName',get(gca,'FontName'),'FontSize',get(gca,'FontSize'));


%%
pltdir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/plots/'];
pltfile = sprintf('%sknm1_NuMassSignal.pdf',pltdir);
export_fig(pltfile);
