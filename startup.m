% this script will be executed every time when Matlab is started
thispath = pwd;
if contains(thispath,'iwsatlas') %when on server
    addpath(genpath('../../../Samak3.0'));
    addpath(genpath('../../Samak3.0'));
    addpath(genpath('../Samak3.0'));
else
    GetSamakPath;
    addpath(genpath(getenv('SamakPath')));
end
