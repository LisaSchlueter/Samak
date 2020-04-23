% statistical neutrino mass sensitivity when background points are excluded
range = 40;
NPfactor = 1.112;
mypath = [getenv('SamakPath'),'knm2ana/knm2_background/results/'];

savename = sprintf('%sknm2_StatSensitivityBackgroundPoints_%.0feV_NPfactor%.2f.mat',mypath,range,NPfactor);

if exist(savename,'file')
    load(savename)
else
    RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
        'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
        'DataType','Twin',...
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2A20',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',NPfactor,...
        'TwinBias_Q',18573.56,...
        'ROIFlag','Default',...
        'chi2','chi2Stat'};
    
    MC = MultiRunAnalysis(RunAnaArg{:});
    MC.exclDataStart = MC.GetexclDataStart(range);
    
    nBkg = 5; % 5 background points
    mNuSq = zeros(nBkg,1);
    mNuSqErr = zeros(nBkg,1);
    
    for i=0:4
        progressbar((i+1)/nBkg);
        MC.exclDataStop = MC.ModelObj.nqU-i;
        MC.Fit;
        
        mNuSq(i+1)    = MC.FitResult.par(1);
        mNuSqErr(i+1) = 0.5*(abs(MC.FitResult.errNeg(1))+MC.FitResult.errPos(1));
    end
    
    qU = MC.RunData.qU;
    
    save(savename,'mNuSq','mNuSqErr','qU','nBkg');
end

%% plot
fbkg = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
x = flipud(qU((end-nBkg+1):end))-18574;
xinterp = linspace(min(x),max(x),1e2);
%plot(xinterp,smooth(interp1(x,mNuSqErr,xinterp,'lin')),'-','LineWidth',2,'Color',rgb('DodgerBlue'));
%hold on;
plot(x,mNuSqErr,'-.o','LineWidth',2,'Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('DodgerBlue'),'MarkerSize',8);
xlabel(sprintf('Upper fit boundary above {\\itE}_0 (eV)'));
ylabel(sprintf('\\sigma({\\itm}_\\nu^2) (eV^2)'));
PrettyFigureFormat('FontSize',22);
t = title(sprintf('NP factor = %.2f',NPfactor),'FontWeight','normal');
xlim([min(x)-5 max(x)+5])
%save
plotdir = strrep(savedir,'results','plots');
MakeDir(plotdir);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.png');
print(plotname,'-dpng','-r500');
fprintf('save plot to %s \n',plotname);