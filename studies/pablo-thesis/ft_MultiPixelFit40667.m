addpath(genpath('../../../Samak2.0'));
close all;
% Select runs to analyze (from 40667 to 40693 are 3h FT runs)
RunList = 40667;
ranges = [9,7,1];
Nranges = length(ranges);

PixelList = [1:100,102:112,114:123];

% Choose fixed parameters (optional)
fixPar = num2str([1,2+2*length(PixelList)+1,2+2*length(PixelList)+2]);

% Choose fitter
fitter = 'matlab';

R40667PreResults = load('data/R40667MultiPix.mat');

for uncertainty = 1:1
    % Choose type of fit (chi2)
    switch uncertainty
        case 1
            chi2name = 'chi2P';
        case 2
            chi2name = 'chi2CM';
    end
    
    for range = 1:Nranges
        % Choose data range to analyze
        dataStart = ranges(range);
        
        switch range
            case 1
                i_N40667 = R40667PreResults.ResultShortStat(3+length(PixelList):3+2*length(PixelList)-1);
                i_B40667 = R40667PreResults.ResultShortStat(3:3+length(PixelList)-1);
            case 2
                i_N40667 = R40667PreResults.ResultMedStat(3+length(PixelList):3+2*length(PixelList)-1);
                i_B40667 = R40667PreResults.ResultMedStat(3:3+length(PixelList)-1);
            case 3
                i_N40667 = R40667PreResults.ResultLongStat(3+length(PixelList):3+2*length(PixelList)-1);
                i_B40667 = R40667PreResults.ResultLongStat(3:3+length(PixelList)-1);
        end
        
        MRA = MultiRunAnalysis('AnaFlag','MultiPixel',...
            'chi2',chi2name,'RunList',RunList,...
            'fitter',fitter,'exclDataStart',dataStart,...
            'fixPar',fixPar,'PixList',PixelList,...
            'DataEffcorr','OFF','i_N',i_N40667,'i_B',i_B40667); %'ROI+PileUp'
        
        MRA.Fit();
        
        mnu = MRA.FitResult.par(1);
        e0 = MRA.FitResult.par(2);
        bck = MRA.FitResult.par(3:3+length(PixelList)-1);
        norm = MRA.FitResult.par(3+length(PixelList):end-2);
        mnu_e = MRA.FitResult.err(1);
        e0_e = MRA.FitResult.err(2);
        bck_e = MRA.FitResult.err(3:3+length(PixelList)-1);
        norm_e = MRA.FitResult.err(3+length(PixelList):end-2);
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
                        ResultShortSys = GeneralResult;
                        
                    case 2
                        ResultMedSys = GeneralResult;
                        
                    case 3
                        ResultLongSys = GenerlaResult;
                end
        end
        
        
        
        
    end
    
    
    
end

save('data/R40667MultiPix');

