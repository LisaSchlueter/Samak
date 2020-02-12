addpath(genpath('../../../Samak2.0'));
close all;

% Select runs to analyze (from 40667 to 40693 are 3h FT runs)
RunList = 'StackCD100all';

%RunList = RunList(1:10);
Nruns = length(RunList);

ranges = [9,7,1];
Nranges = length(ranges);
pixellist = 1;

% Apply any pulls you wish (recommended for FT: 2 on mnu, none on the others)
pulls    = [Inf,Inf,Inf,Inf,Inf,Inf];

% Choose fixed parameters (optional)
fixPar   = '1 5 6';

% Choose fitter
fitter = 'minuit';

for uncertainty = 1:2
    % Choose type of fit (chi2)
    switch uncertainty
        case 1
            chi2name = 'chi2Stat';
        case 2
            chi2name = 'chi2CM';
    end
    
    for range = 1:Nranges
        % Choose data range to analyze
        dataStart = ranges(range);
        
        GeneralResult = zeros(Nruns,10);
        MRA = MultiRunAnalysis('AnaFlag','StackPixel',...
            'chi2',chi2name,'RunList',RunList,... % RunList, 
            'fitter',fitter,'exclDataStart',dataStart,...
            'fixPar',fixPar,'pulls',pulls,...
            'DataEffcorr','OFF',...
            'EffCorrFlag','ON'); %'ROI+PileUp'
        
        MRA.Fit();
        
        mnu = MRA.FitResult.par(1);
        e0 = MRA.FitResult.par(2);
        bck = MRA.FitResult.par(3:3+length(pixellist)-1);
        norm = MRA.FitResult.par(3+length(pixellist):end-2);
        mnu_e = MRA.FitResult.err(1);
        e0_e = MRA.FitResult.err(2);
        bck_e = MRA.FitResult.err(3:3+length(pixellist)-1);
        norm_e = MRA.FitResult.err(3+length(pixellist):end-2);
        chi2 = MRA.FitResult.chi2min;
        dof = MRA.FitResult.dof;
        GeneralResult(1,:) = [mnu,e0,bck,norm,mnu_e,e0_e,bck_e,norm_e,chi2,dof];
        
        
        
        switch uncertainty
            case 1 % Statistical
                switch range
                    case 1
                        ResultShortStat = GeneralResult;
                        
                    case 2
                        ResultMedStat = GeneralResult;
                        
                    case 3
                        ResultLongStat = GeneralResult;
                        
                end
            case 2 % systematics with CovMat
                switch range
                    case 1
                        fprintf('-----Short range above----\n');
                        ResultShortSys = GeneralResult;
                        
                    case 2
                        fprintf('-----Medium range above----\n');
                        ResultMedSys = GeneralResult;
                        
                    case 3
                        fprintf('-----Long range above-----\n');
                        ResultLongSys = GeneralResult;
                end
        end
 
    end

end

save('data/StackedPixelsStackedRunsResultsThesis2');






