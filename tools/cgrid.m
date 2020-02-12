%------------------------------------------------------------------------------------
% PROGRAM NAME: cgrid
% AUTHOR: Andri m. Gretarsson
% DATE: Jun 3, 2004
%
% SYNTAX: <gridaxeshandle>=cgrid(<colorarg, stylearg, colorargy>)
%           <...> indicates optional arguments
%
% Creates a new set of axes containing the colored grid, 
% immediately behind the current axes.  To make the axes with the
% grid visible, the current axes are made transparent and the
% current axes grid is removed.  If any colored grid axes are
% in use when this function is called, it first deletes these
% old grids, before adding new cgrid axes.  That way figure children
% don't multiply.  
% 
%
% INPUT ARGUMENTS:
% ----------------
% colorarg  = The color of the grid.  The form of colorarg may be
%             any of the following forms:
%
%             'g'               (or any single-letter color)
%             [0.3 0.5 0.6]     (or any rgb combination)
%             'g--'             (or any color-style combination having the
%                                same syntax as in the plot matlab function.
%
% stylearg  = The linestyle of the grid lines.  Format is one of:
%             '-',':' or '--'   (supercedes any style specified in colorarg)
%
% colorargy = The color of the vertical grid lines. Horizontal grid lines
%             will remain the color specified in colorarg.  This argument has
%             the same format as colorarg except that it may not include the
%             linestyle specification (for example, 'g' or [0 1 0] is OK but 
%             not 'g-').
%
%
% OUTPUT ARGUMENTS:
% -----------------
% gridhandle = The handle to the set of axes containing the colored
%              grid.  This can be used to set the background color
%              and many other things.
%
%
% NOTE:
% -----
%
% Although this function correctly treats legends already present when the 
% function is called, legends that are added later will have the same
% color background as the figure has when the legend is called (normally grey).
% This is due to the necessity of making the current figure transparent so 
% that the cgrid is visible. However, this effect does NOT show up on exported 
% plots.  If you are unhappy with the transparent legend bacground in the
% matlab plot itself, this can be changed by saving the legend handle as
% in the following example and setting the color to the desired color.
%
% EXAMPLE:
% --------
%       legendhandle=legend('test legend');
%       set(legendhandle,'color','w');
%
%
% KNOWN ISSUES:
% -------------
% Unlike the built-in grid function, The number of cgrid lines does not change when
% the figure window size is changed.  If this is an issue, users should set the figure
% window size before applying cgrid.
%
% LAST UPDATED: 06/28/04 by AMG.
%------------------------------------------------------------------------------------
% SYNTAX: <gridaxishandle>=cgrid(<colorarg, stylearg, colorargy>)
%------------------------------------------------------------------------------------

function varargout=cgrid(varargin);

defaultcolor=[0.8 0.8 0.8];
defaultstyle=':';

if nargin==0
    theXcolor=defaultcolor;
    thestyle=defaultstyle;
end

if nargin>=1                                    % Parse the input arguments
    colorarg=varargin{1};
    if length(colorarg)==1
        if sum(colorarg(1)=='brgmcykw')>=1 
            theXcolor=colorarg(1);
            thestyle=defaultstyle;
        else
            theXcolor=[0.8 0.8 0.8];
            thestyle=colorarg(1);
        end
    else
        if length(colorarg)==2
            if sum(colorarg(1)=='brgmcykw')>=1
                theXcolor=colorarg(1);
                thestyle=colorarg(2);
            else
                theXcolor=defaultcolor;
                thestyle=colorarg;
            end
        else
            if ischar(colorarg(3))
                theXcolor=colorarg(1);
                thestyle=colorarg(2:3);
            else
                theXcolor=colorarg(1:3);
                thestyle=defaultstyle;
            end
        end
    end
end
theYcolor=theXcolor;

if nargin>=2 thestyle=varargin{2}; else style=':'; end
if nargin>=3 theYcolor=varargin{3}; end

grid off;
oldaxes=gca;
oldaxesprop=get(gca);
set(oldaxes,'color','none');                    % remove white background
set(oldaxes,'xgrid','off');                     % remove grid
set(oldaxes,'ygrid','off');                     % remove grid
set(oldaxes,'xminorgrid','off');                % remove grid
set(oldaxes,'yminorgrid','off');                % remove grid

kids=get(gcf,'children');                       % remove any previous cgrid axes
s=1;
while s<=length(kids)
    tag=get(kids(s),'tag');
    tagstring=['gridonly',num2str(oldaxes,'%0.4f')];
    if length(tag)==length(tagstring) & tag==tagstring;
        delete(kids(s));
        disp('Removed old cgrid');
    end
    s=s+1;
end

newaxes=axes;                                   % create axes for the colored grid
set(newaxes,'Position',oldaxesprop.Position);
kids=get(gcf,'children');                       % read all the figure axes
oldkid=find(kids==oldaxes);                     % insert the new cgrid axes           
kids=[kids(2:oldkid);newaxes;kids(oldkid+1:end)]; % behind the set of axes that was 
set(gcf,'children',kids);                       % current when the function was called
                                                
set(newaxes,'xlim',get(oldaxes,'xlim'));
set(newaxes,'ylim',get(oldaxes,'ylim'));
set(newaxes,'xtick',get(oldaxes,'xtick'));
set(newaxes,'ytick',get(oldaxes,'ytick'));
set(newaxes,'xscale',get(oldaxes,'xscale'));
set(newaxes,'yscale',get(oldaxes,'yscale'));
set(newaxes,'xgrid','on');
set(newaxes,'ygrid','on');
set(newaxes,'gridlinestyle',thestyle);
set(newaxes,'xcolor',theXcolor);
set(newaxes,'ycolor',theYcolor);
set(newaxes,'tag',...                           % apply tag so cgrid axes can be recognized. Need the
    ['gridonly',num2str(oldaxes,'%0.4f')]);     % oldaxes handle so that cgrids on other subplots don't get erased.
                                                
set(newaxes,'xticklabel',[]);
set(newaxes,'yticklabel',[]);
set(newaxes,'box','off');

if nargout>=1, varargout{1}=newaxes; end
set(gcf,'currentaxes',oldaxes);

