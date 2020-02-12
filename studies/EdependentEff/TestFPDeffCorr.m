% Option: data/sim
Mode     = 'Data';
%Mode     = 'Sim';

% Run List
RunList = 0;
if RunList == 0
    [ n, datafiles , RunList ] = GetRunList( '../../tritium-data/hdf5/','*.h5',1,'string');
    RunList=RunList(RunList>40666 & RunList<40694);
end

% Run Analysis
chi2='chi2Stat';
A = MultiRunAnalysis('RunList',RunList,'chi2',chi2,'exclDataStart',1,...
    'fixPar','1 5 6','ringCutFlag','ex2','AnaFlag','StackPixel');

% Option for simulation
switch Mode
    case 'Sim'
        A.RunData.TBDIS = zeros(A.ModelObj.nqU,1);
        
        A.SingleRunObj = cell(length(A.StackedRuns),1);
        for r=1:length(A.StackedRuns)
            A.SingleRunObj{r} = ref_RunSummaries_StackPix(A.StackedRuns(r),A.ringCutFlag,...
                'ISCS','Theory','recomputeRF','OFF');
            A.SingleRunObj{r}.ComputeTBDDS;
            A.SingleRunObj{r}.ComputeTBDIS;
        end
        
        for i=1:1:numel(RunList)
            %A.SingleRunObj{i}.AddStatFluctTBDIS;
            A.RunData.TBDIS  = A.RunData.TBDIS + A.SingleRunObj{i}.TBDIS;
        end
        A.RunData.qU     = A.ModelObj.qU;
end

% Corrections
% Correct Data for ROI Efficiency
ROIEffCorrectionFactor = @(qu) -1.139949799271709e-04*(qu/1e3).^4 + 0.00773350959014*(qu/1e3).^3 + -0.196747862772578*(qu/1e3).^2 + 2.227682510553492*(qu/1e3) -8.488123971598990;            
ROI  = 1./ROIEffCorrectionFactor(A.RunData.qU);
% Correct Data for PileUp Efficiency
PUeffCorrectionFactor = @(pixelrate) -1.597653658180098e-06.*pixelrate + 0.999999590625725;
PU   = 1./PUeffCorrectionFactor(A.RunData.TBDIS./A.RunData.qUfrac./A.RunData.TimeSec/124);
%% Plot
figure(1);
plot(A.RunData.qU,A.RunData.TBDIS./A.RunData.TBDIS,A.RunData.qU,ROI,A.RunData.qU,PU,A.RunData.qU,ROI.*PU);
plt = Plot();
plt.XLabel = 'qU (V)'; % xlabel
plt.YLabel = 'Multiplicative Correction'; %ylabel
plt.Title = 'ROI and Pile-Up Efficiency Corrections to FT Data'; % plot title
plt.YGrid = 'on'; 
plt.XGrid = 'on'; 
plt.YScale = 'linear';
plt.XScale = 'linear';
plt.Legend =  {'No correction', 'ROI', 'Pile-Up','ROI+Pile-Up'}; 
plt.BoxDim = [7, 3]; %[width, height] in inches
plt.export('ROI-PileUp-EfficencyCorrection.png');

%% Apply Correction or Not
A.RunData.TBDIS = A.RunData.TBDIS.*PU.*ROI;

%% Fit
A.Fit('CATS','OFF');
A.PlotFit;