addpath(genpath('../../../Samak2.0'));
close all;
    
    RunList = 'StackCD100all';
    ranges = [9,7];
    Nranges = length(ranges);
    
    PixelList = [1:100,102:112,114:123];
    
    % Choose fixed parameters (optional)
    fixPar = num2str([1,2+2*length(PixelList)+1,2+2*length(PixelList)+2]);
    
    % Choose fitter
    fitter = 'matlab';

%Mpix = load('data/MultiPixelResultsThesis2.mat');

for uncertainty = 2:2
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
        
        switch range
            case 11
                MpixN = Mpix.ResultShortStat(3+length(PixelList):3+2*length(PixelList)-1);
                MpixB = Mpix.ResultShortStat(3:3+length(PixelList)-1);
                MpixQ = Mpix.ResultShortStat(2);
            case 21
                MpixN = Mpix.ResultMedStat(3+length(PixelList):3+2*length(PixelList)-1);
                MpixB = Mpix.ResultMedStat(3:3+length(PixelList)-1);
                MpixQ = Mpix.ResultMedStat(2);
%                   MpixN = []; MpixB = []; MpixQ = 0;

        end
        MpixN = []; MpixB = []; MpixQ = 0;
        MRA = MultiRunAnalysis('AnaFlag','MultiPixel',...
            'chi2',chi2name,'RunList',RunList,...
            'fitter',fitter,'exclDataStart',dataStart,...
            'fixPar',fixPar,'PixList',PixelList,...
            'DataEffcorr','OFF','EffCorrFlag','ON',...
            'i_N',MpixN,'i_B',MpixB,'i_Q',MpixQ); %'ROI+PileUp'
        
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
clear Mpix;
save('data/MultiPixelResultsThesis2');

