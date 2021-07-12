
function [AsymErrNeg, AsymErrPos] = GetAsymmErrMinos(filename)
% ---------------------------------------------------------------------------------
% get asymmetric uncertainties from MINOS output
% inside fit class (FITC), minuit command window output is stored in text file 
% this script:
% - reads this text file
% - searches for pattern for neutrino mass, endpoint, Bkg, N fit results
% output: asymmetric minos uncertainties
%
% L. Schlüter, Feb. 2020
% ---------------------------------------------------------------------------------
d = fileread(filename);
minuitStr = 'ERROR      NEGATIVE      POSITIVE';
ResultIndex = strfind(d,minuitStr);
dshort = d(ResultIndex(end):end);

% look for patters
FindmNuStr = ' 1   Par  #1'; % neutrino mass
FindE0Str  = ' 2   Par  #2'; % endpoint
FindBkgStr = ' 3   Par  #3'; % background
FindNStr   = ' 4   Par  #4'; % uniform
FindEndStr = ' 5   Par  #5'; % end str
FindetaStr = '18   Par #18'; % relic neutrino overdensity

% find index in text
mNuStartIndex = strfind(dshort,FindmNuStr);
E0StartIndex  = strfind(dshort,FindE0Str);
BkgStartIndex = strfind(dshort,FindBkgStr);
NStartIndex   = strfind(dshort,FindNStr);
EndIndex      = strfind(dshort,FindEndStr);
etaStartIndex = strfind(dshort,FindetaStr);

%% extract asymmetric uncertainties
AsymErrNeg = zeros(5,1);
AsymErrPos = zeros(5,1);

try
    [AsymErrNeg(1), AsymErrPos(1)] = GetErr(dshort,mNuStartIndex,E0StartIndex,FindmNuStr);
catch
end

try
catch
    [AsymErrNeg(2), AsymErrPos(2)] = GetErr(dshort,E0StartIndex,BkgStartIndex,FindE0Str);
end

try
    [AsymErrNeg(3), AsymErrPos(3)] = GetErr(dshort,BkgStartIndex,NStartIndex,FindBkgStr);
catch
end

try
    [AsymErrNeg(4), AsymErrPos(4)] = GetErr(dshort,NStartIndex,EndIndex,FindNStr);
catch
end
try
    [AsymErrNeg(5), AsymErrPos(5)] = GetErr(dshort,etaStartIndex,etaStartIndex+70,FindetaStr);
catch
end

    function [errNeg, errPos] = GetErr(dshort,StartIndex,StopIndex,FindStr)
        Mystr = extractAfter(dshort(StartIndex:StopIndex),FindStr);
        if contains(Mystr,'fixed')
            errNeg = 0;
            errPos = 0;
        else
            FitResultNum = str2num(Mystr);
            errNeg = FitResultNum(3);
            errPos = FitResultNum(4);
        end
    end
end
