addpath(genpath('../../../Samak2.0'));
close all;


    % Select runs to analyze (from 40667 to 40693 are 3h FT runs)
    
    RunList = load('RunList100Good.mat');
    RunList = RunList.RunList100Good;
    RunList((RunList == 40610) | (RunList == 40612) | (RunList == 40539)) = [];
    

    %RunList = RunList(1:10);
    Nruns = length(RunList);
    
    ranges = [9,7,1];
    Nranges = length(ranges);pixellist = 1;
    
    % Apply any pulls you wish (recommended for FT: 2 on mnu, none on the others)
    pulls    = [Inf,Inf,Inf,Inf,Inf,Inf];
    
    % Choose fixed parameters (optional)
    fixPar   = '1 5 6';
    
    % Choose fitter
    fitter = 'minuit';


for uncertainty = 1:1
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
        
        GeneralTable = zeros(Nruns,10);
        parfor r = 1:Nruns
            run = RunList(r);
            
            RA = RunAnalysis('AnaFlag','StackPixel',...
                'chi2',chi2name,'RunNr',run,...
                'fitter',fitter,'exclDataStart',dataStart,...
                'fixPar',fixPar,'pulls',pulls,...
                'DataEffcorr','OFF'); %'ROI+PileUp'
            
            RA.Fit();
            fprintf('Run: %d (%g / %g) \n',run,r,Nruns);
            
            mnu = RA.FitResult.par(1);
            e0 = RA.FitResult.par(2);
            bck = RA.FitResult.par(3:3+length(pixellist)-1);
            norm = RA.FitResult.par(3+length(pixellist):end-2);
            mnu_e = RA.FitResult.err(1);
            e0_e = RA.FitResult.err(2);
            bck_e = RA.FitResult.err(3:3+length(pixellist)-1);
            norm_e = RA.FitResult.err(3+length(pixellist):end-2);
            chi2 = RA.FitResult.chi2min;
            dof = RA.FitResult.dof;
            GeneralTable(r,:) = [mnu,e0,bck,norm,mnu_e,e0_e,bck_e,norm_e,chi2,dof];
            
            
        end
        
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
        
        
        
        
    end
    
    
    
end


save('StackedPixelsResultsThesis');