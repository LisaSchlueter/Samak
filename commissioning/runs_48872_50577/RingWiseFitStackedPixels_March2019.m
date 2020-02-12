%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ring fit analysis
% ---------------------------------------------------------------------- %
% First Tritium Data ring-wise, for stacked runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../../../Samak2.0'));

    % Select runs to analyze 
    RunList = 'March2019';
    % Select the rings you want to analyze, ex. [1,5,8], [4:10], etc.
    ringlist = [1:12];

    % Apply any pulls you wish 
    pulls    = [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf];
    
    % Choose data range to analyze
    %ranges = [2 8 14 20];  
    ranges = [2 14];  
    Nranges = length(ranges);
    
    % Choose fixed parameters (optional)
    fixPar   = '1 5 6 7 8 9 10';
    
    % Choose fitter
    fitter = 'minuit';

for uncertainty = 1
    % Choose type of fit (chi2)
    switch uncertainty
        case 1
            chi2name = 'chi2Stat';
        case 2
            chi2name = 'chi2CMShape';
    end
    
    for range = 1:Nranges
        
        dataStart = ranges(range);
        
        % Call Analysis class
        MRA = MultiRunAnalysis(...
            'DataType','Real',...
            'ELossFlag','KatrinD2',...
            'StackTolerance',1,...
            'NonPoissonScaleFactor',1.5,...
            'AnaFlag','Ring','RingList',ringlist,...
            'chi2',chi2name,'RunList',RunList,...
            'fitter',fitter,'exclDataStart',dataStart,...
            'fixPar',fixPar,'pulls',pulls);
        
        % Do the fits
        MRA.FitAllRings();
        
        % Save the results for the plotting class
        FitResult.par     = MRA.AllFitResults{1};
        FitResult.err     = MRA.AllFitResults{2};
        FitResult.chi2min = MRA.AllFitResults{3};
        FitResult.dof     = MRA.AllFitResults{5};
        
        GeneralTable = [FitResult.par,FitResult.err,FitResult.chi2min,FitResult.dof];
        
        switch uncertainty
            case 1 % Statistical
                switch range
                    case 1
                        TableShortStat = GeneralTable;
                        
                    case 2
                        TableMedStat = GeneralTable;
                        
                    case 3
                        TableLongStat = GeneralTable;
                        
                end
            case 2 % systematics with CovMat
                switch range
                    case 1
                        TableShortSys = GeneralTable;
                        
                    case 2
                        TableMedSys = GeneralTable;
                        
                    case 3
                        TableLongSys = GenerlaTable;
                end
        end
        
        PlotRingWiseFit(MRA,dataStart,GeneralTable);

    end
    
end

%save('results/StackedRings');

%% Plot
function PlotRingWiseFit(obj,dataStart,GeneralTable)

index = obj.ModelObj.Q_i-obj.ModelObj.qU;
myMainTitle=[sprintf('KATRIN Ring-wise Fit - %d Runs (%s) - %.1f',...
    numel(obj.RunList),obj.chi2,index(dataStart)),'eV below E0'];
maintitle=myMainTitle;
savefile=sprintf('plots/runwise/MRA_RunWiseFitResults%d_%s_%.0feVbelowE0-RingWise.pdf',numel(obj.RunList),obj.chi2,index(dataStart));
fig1 = figure('Name','MRA Ring-wise Fits','NumberTitle','off','rend','painters','pos',[10 10 1000 600]);
a=annotation('textbox', [0 0.91 1 0.1], ...
    'String', maintitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

bar(obj.RingList,GeneralTable(:,2)+obj.ModelObj.Q_i-mean(GeneralTable(:,2)+obj.ModelObj.Q_i),'facecolor',rgb('IndianRed'));
hold on
errorb(obj.RingList,GeneralTable(:,2)+obj.ModelObj.Q_i-mean(GeneralTable(:,2)+obj.ModelObj.Q_i),GeneralTable(:,12))
hold off

ylabel(sprintf('E_0-%.1f (eV)',obj.ModelObj.Q_i));
xlabel('Ring');
sp1 = sprintf('<E0>=%.2f eV \\pm %.2f eV (std)',mean(GeneralTable(:,2)+obj.ModelObj.Q_i),std(GeneralTable(:,2)));title(sp1)
grid on
PrettyFigureFormat
a.FontSize=20;a.FontWeight='bold';
publish_figurePDF(gcf,savefile);
end

