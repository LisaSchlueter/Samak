% this script will be executed every time when Matlab is started
thispath = pwd;

if contains(thispath,'iwsatlas') % when on server
% <<<<<<< HEAD
% addpath(genpath('../../../Samak3.0'));
% addpath(genpath('../../Samak3.0'));
% addpath(genpath('../Samak3.0'));
% else
% 
% if contains(thispath,'Samak') % if working with Samak
% GetSamakPath;
% addpath(genpath(getenv('SamakPath')));
% elseif contains(thispath,'roplab') % if working for roplab
% % roplab
% GetRoplabPath;
% addpath(genpath(getenv('RoplabPath')));
% end 
% end
% =======
    addpath(genpath('../../../Samak3.0'));
    addpath(genpath('../../Samak3.0'));
    addpath(genpath('../Samak3.0'));
else
  
    if contains(thispath,'Samak')   % if working with Samak
        GetSamakPath;
        addpath(genpath(getenv('SamakPath')));
    elseif contains(thispath,'roplab')   % if working for  roplab
        % roplab
        GetRoplabPath;
        addpath(genpath(getenv('RoplabPath')));
    end    
end
% >>>>>>> d8ad3705616652c7335ddc80e7a9e273001abd2c
