% Return and edit a colormap based on the predefinde matlab-colormaps. 
% The color map can be rescaled and the number of steps can be adjusted.
% 
% Initialize the colormap
% EditColormap('InitialMap','default');
% 
% Change the colorrange
% EditColormap('Scale',[0 1]);
% 
% Set the number of colors int the gradient.
% EditColormap('Size',10);


function [cmp] = EditColormap(varargin)

p = inputParser;    p.CaseSensitive     = 0;  
                    p.FunctionName      = 'EditColormap.m';
                    p.KeepUnmatched     = 0;
                    p.PartialMatching   = 1;
                    p.StructExpand      = 1;

addParameter(p,'Axis',gca);                    
addParameter(p,'InitialMap','default',@(x)ismember(x,...
    {'default','parula','jet','hsv','hot','cool','spring','summer','autumn','winter',...
    'gray','bone','copper','pink','lines','colorcube','prism','flag','white'}));
addParameter(p,'Scale',[0 1],@(x)assert(isfloat(x) && (~any(x<0)) && (~any(x>1)),...
    'Check if [x y] are in correct limits: 0<=x(:)<=1'));
addParameter(p,'Size',0,@(x)isfloat(x) && x>0);

parse(p,varargin{:});

% use the parameter: InitialMap
colormap(p.Results.InitialMap); cmpt=colormap; 
sv=size(cmpt); s=sv(1);

% use the parameter: Scale
sc1 = p.Results.Scale(1);   sc2 = p.Results.Scale(2);

if sc1 < sc2
    left = 1+sc1*s;   right = sc2*s;
elseif sc1 > sc2
    left = sc1*s;     right = 1+sc2*s;
elseif sc1 == 1
    left = s;         right = s;
elseif sc1 ==0
    left = 1;            right = 1;
else
    left = .5+sc1*s;  right = .5+sc1*s;
end

% use the paramter: Size
st = s;
if p.Results.Size > 0
    st = int8(p.Results.Size);
end

% calculate final colormap
cmpt2=interp1(1:s(1),cmpt, linspace(left,right,st) ); 

% change the colormap
cmp=cmpt2; 

clear cmpt2 sc1 sc2 left right st s;
colormap(p.Results.Axis,cmp);
