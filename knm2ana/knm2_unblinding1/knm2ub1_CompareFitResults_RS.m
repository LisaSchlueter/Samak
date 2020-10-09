
path = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/results/'];
name = [path,'knm2ub1_CompareRSfitresults.mat'];

if exist(name,'file') 
    load(name,'E0','N','B','chi2min','E0err','Berr','Nerr','LiveTime');
else
    
    range = 40;
    RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
        'fixPar','E0 Bkg Norm',...           % free Parameter !!
        'DataType','Real',...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...
        'chi2','chi2Stat'};
    
    A = MultiRunAnalysis(RunAnaArg{:});
    A.exclDataStart = A.GetexclDataStart(range);
    %%
    path1 = [getenv('SamakPath'),'tritium-data/fit/Knm2/Uniform/'];
    path2 = [getenv('SamakPath'),'tritium-data/fit/Knm2/Uniform_RS5g/'];
    
    E0       = zeros(A.nRuns,2);   E0err    = zeros(A.nRuns,2);
    B        = zeros(A.nRuns,2);   Berr        = zeros(A.nRuns,2);
    N        =  zeros(A.nRuns,2);  Nerr        =  zeros(A.nRuns,2);
    chi2min  =  zeros(A.nRuns,2);
    
    for i=1:A.nRuns
        progressbar(i/A.nRuns);
        nameRS_5i = sprintf('%sFit%.0f_chi2Stat_39bE0_freeParE0BkgNorm.mat',path1,A.RunList(i));
        d1           =  importdata(nameRS_5i);
        E0(i,1)      = d1.FitResult.par(2)+A.ModelObj.Q_i;
        N(i,1)       = d1.FitResult.par(4)+1;
        B(i,1)       = d1.FitResult.par(3)+A.ModelObj.BKG_RateSec_i;
        chi2min(i,1) = d1.FitResult.chi2min;
        E0err(i,1)   = d1.FitResult.err(2);
        Nerr(i,1)    = d1.FitResult.err(4);
        Berr(i,1)    = d1.FitResult.err(3);
        
        nameRS_5g = sprintf('%sFit%.0f_chi2Stat_40bE0_freeParE0BkgNorm.mat',path2,A.RunList(i));
        d2           = importdata(nameRS_5g);
        E0(i,2)      = d2.FitResult.par(2)+A.ModelObj.Q_i;
        N(i,2)       = d2.FitResult.par(4)+1;
        B(i,2)       = d2.FitResult.par(3)+A.ModelObj.BKG_RateSec_i;
        chi2min(i,2) = d2.FitResult.chi2min;
        E0err(i,2)   = d2.FitResult.err(2);
        Nerr(i,2)    = d2.FitResult.err(4);
        Berr(i,2)    = d2.FitResult.err(3);
    end
    
    RSversion = {'Durable-5i','Durable-5g'};
    
    LiveTime = hours(A.SingleRunData.StartTimeStamp-A.SingleRunData.StartTimeStamp(1));
     
    MakeDir(path)
    save(name,'E0','N','B','chi2min','RSversion','RunAnaArg','LiveTime');
end

%% plot
 PlotStyle = { '.','MarkerSize',25};
 plotpath = [getenv('SamakPath'),'knm2ana/knm2_unblinding1/plots/'];
MakeDir(plotpath);
 %% E0 
fig88 = figure('units','normalized','pos',[0.1, 0.1,1.0,0.5]);
pline = plot(linspace(-5,1500,100),mean(E0(:,1)-E0(:,2)).*ones(100,1),'k:','LineWidth',2); 
hold on;
e1 = plot(LiveTime,E0(:,1)-E0(:,2),PlotStyle{:},'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',22);
xlim([-0.5,max(LiveTime)+0.5])
ylabel(sprintf('{\\itE}_0^{RS-5i} - {\\itE}_0^{RS-5g} (eV)'));
xlabel('Time (hours)');
leg = legend([e1,pline],'KNM2 scan-wise fits',sprintf('\\mu = %.3f eV , \\sigma = %.3f eV',mean(E0(:,1)-E0(:,2)),std(E0(:,1)-E0(:,2))),'EdgeColor',rgb('Silver'),'Location','southeast');
t = title('Endpoint','FontWeight','normal','FontSize',get(gca,'FontSize'));
plotname1 = [plotpath,'knm2ub1_CompareRSfitresults_E0.png'];
print(plotname1,'-dpng','-r350');
%% background
fig88 = figure('units','normalized','pos',[0.1, 0.1,1.0,0.5]);
pline = plot(linspace(-5,1500,100),mean((B(:,1)-B(:,2))*1e3).*ones(100,1),'k:','LineWidth',2); 
hold on;
e1 = plot(LiveTime,(B(:,1)-B(:,2))*1e3,PlotStyle{:},'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',22);
xlim([-0.5,max(LiveTime)+0.5])
ylabel(sprintf('{\\itB}^{RS-5i} - {\\itB}^{RS-5g} (mcps)'));
xlabel('Time (hours)');
leg = legend([e1,pline],'KNM2 scan-wise fits',sprintf('\\mu = %.1g mcps , \\sigma = %.1g mcps',...
    mean((B(:,1)-B(:,2))*1e3),std((B(:,1)-B(:,2))*1e3)),'EdgeColor',rgb('Silver'),'Location','southeast');
t = title('Background','FontWeight','normal','FontSize',get(gca,'FontSize'));
plotname2 = [plotpath,'knm2ub1_CompareRSfitresults_B.png'];
print(plotname2,'-dpng','-r350');
%% normalization
fig88 = figure('units','normalized','pos',[0.1, 0.1,1.0,0.5]);
pline = plot(linspace(-5,1500,100),mean(N(:,1)-N(:,2)).*ones(100,1),'k:','LineWidth',2); 
hold on;
e1 = plot(LiveTime,N(:,1)-N(:,2),PlotStyle{:},'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',22);
xlim([-0.5,max(LiveTime)+0.5])
ylabel(sprintf('{\\itN}^{RS-5i} - {\\itN}^{RS-5g}'));
xlabel('Time (hours)');
leg = legend([e1,pline],'KNM2 scan-wise fits',sprintf('\\mu = %.1g , \\sigma = %.1g',...
    mean(N(:,1)-N(:,2)),std(N(:,1)-N(:,2))),'EdgeColor',rgb('Silver'),'Location','southeast');
t = title('Signal normalization','FontWeight','normal','FontSize',get(gca,'FontSize'));
plotname3 = [plotpath,'knm2ub1_CompareRSfitresults_N.png'];
print(plotname3,'-dpng','-r350');

%% chi2min
fig88 = figure('units','normalized','pos',[0.1, 0.1,1.0,0.5]);
dof = 25;
pVal = 1-chi2cdf(chi2min,dof);
pline = plot(linspace(-5,1500,100),mean(pVal(:,1)-pVal(:,2)).*ones(100,1),'k:','LineWidth',2); 
hold on;
e1 = plot(LiveTime,pVal(:,1)-pVal(:,2),PlotStyle{:},'Color',rgb('DodgerBlue'));
PrettyFigureFormat('FontSize',22);
xlim([-0.5,max(LiveTime)+0.5])
ylabel(sprintf('{\\itp}^{RS-5i} - {\\itp}^{RS-5g}'));
xlabel('Time (hours)');
leg = legend([e1,pline],'KNM2 scan-wise fits',sprintf('\\mu = %.1g , \\sigma = %.1g',...
    mean(pVal(:,1)-pVal(:,2)),std(pVal(:,1)-pVal(:,2))),'EdgeColor',rgb('Silver'),'Location','southeast');
t = title('p-value','FontWeight','normal','FontSize',get(gca,'FontSize'));
plotname3 = [plotpath,'knm2ub1_CompareRSfitresults_pVal.png'];
print(plotname3,'-dpng','-r350');