addpath(genpath('../../../Samak2.0'));

RunList = 'StackCD100all';

ranges = [9,7,1]; % 9 = -200 eV; 7 = -400 eV; 1 = -1600 eV
Nranges = length(ranges);

% Choose fixed parameters (optional)
fixPar   = '1 5 6';

% Choose fitter
fitter = 'minuit';

% Choose CM 
chi2name = 'chi2CM';


for range = 1:Nranges
    % Choose data range to analyze
    dataStart = ranges(range);
    
    MRA = MultiRunAnalysis('AnaFlag','StackPixel',...
        'chi2',chi2name,'RunList',RunList,... 
        'fitter',fitter,'exclDataStart',dataStart,...
        'fixPar',fixPar,...
        'DataEffcorr','OFF');
    
    MRA.Fit();
    
end



