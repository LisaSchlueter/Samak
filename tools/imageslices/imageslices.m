function [axis_handles, yslices] = imageslices(varargin)
% This function plots matrix data using imagesc, with subplots showing
% horizontal slices through the data. Each call to imageslices creates a
% new MATLAB figure.
%
% Syntax:
% imageslices(z)
% imageslices(x, y, z)
% imageslices(z, Nslices)
% imageslices(z, Nslices, ylimits)
% imageslices(z, ylimits)
% imageslices(___, Name, Value)
% imageslices(___, Property)
%
% Inputs:
%   z: 2D matrix of data values, size M-by-N (required)
%   x: Vector of x values, must have length N (optional)
%   y: Vector of y values, must have length M (optional)
%   Nslices: Number of slices to take (optional, default = 10)
%   ylimits: 2-element vector with coordinates of the lowest & highest
%               slices. If a y vector is provided, these coordinates will  
%               be rounded to the nearest values present in y. If a y
%               vector is not provided, these coordinates represent row
%               indexes in z. (optional, default = [min(y),max(y)]).
%   Name, Value: Specifies additional plotting properties using one or more
%                   name-value pair arguments. These are described below.
%   Property: Specifies additional plotting properties that are not in
%                   name-value pair arguments.
%
% Outputs:
%   1) The newly created figure
%   2) axis_handles: a vector of handles for the various axes in the new
%         figure, starting with the imagesc axis.
%   3) yslices: vector of row indexes at which the slices of z were taken
%
% Name-Value Property Pairs:
%   'labels': Cell array of labels for the X,Y,Z axes. Default = {'x','y','z'}
%   'colors': Colors to use for the slices. Must be a matrix of size
%                Nslices-by-3. Default = lines(Nslices)
%   'ticklength': Length of the colored ticks that are added on top of the
%                    imagesc axis, in axis-normalized units. Default = 0.02
%   'tickstyle': Style of the colored ticks that are added on top of the
%                    imagesc axis; any valid MATLAB line style is allowed.
%                    Default = '-'
%   'topmargin': Margin above the top-most subplot of data slices, in
%                   figure-normalized units. Default = 0.05
%   'bottommargin': Margin below the bottom-most subplot of data slices, in
%                     figure-normalized units. Default = 0.1
%   'subplotgap': Vertical gap between the subplots of data slices, in
%                    figure-normalized units. Default = 0.005
%   'deltaz': Including this property causes all slices to be plotted in a
%                single subplot. The numeric Value following this property
%                Name will be used as a vertical spacing between data slices.
%
% Additional Plotting Properties:
%   'subtractmean': This property will only take effect if 'deltaz' is also
%                      in use. Before each data slice is plotted, its mean
%                      will be subtracted.
%
% Notes:
%   If slices are to be plotted in separate subplots, each subplot
%      will have the save vertical axis limits; these will by default be
%      min(z) and max(z).
%   If more than 3 subplots are created for the data slices, yticklabels
%      will only be shown for every other subplot.
%
%
% Examples:
% imageslices(peaks(100));
% imageslices(peaks(100), 5, [20,80]);
% imageslices(peaks(100), [20,80], 'deltaz', 5);
% imageslices(peaks(100), [20,80], 'deltaz', 5, 'subtractmean');
% imageslices(peaks(100), 5, [20,80], 'ticklength', 1, 'tickstyle', '--');
% x = (1:100);
% y = (1:100);
% z = peaks(100);
% hax = imageslices(x, y, z, 'labels', {'x label','y label','z label'});
% colormap(hax(1), 'jet');
%
%
% v1.2, Bob De Alba
% 2018-03-20: original release
% 2018-03-20: Corrected the vertical positioning of colored ticks on the
%             imagesc axis to account for the imagesc pixel size.
% 2018-03-21: Changed the format and ordering of input variables. Added the
%             option to show all data slices in a single subplot via the
%             'deltaz' input. Also added the 'ticklength', 'tickstyle',
%             'topmargin', 'bottommargin', 'subplotgap', and 'subtractmean'
%             properties. Added the 'yslices' output. Rewrote the help
%             documentation.

if all(size(varargin{1}) > 1) % check if 1st input is a vector or matrix
    z = varargin{1};
    x = (1:size(z,2));
    y = (1:size(z,1));
    varargin(1) = [];
else
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    varargin(1:3) = [];
end

% check if the user has specified Nslices
if ~isempty(varargin) && isa(varargin{1},'numeric') && length(varargin{1})==1
    Nslices = varargin{1};
    varargin(1) = [];
else
    Nslices = 10;
end

% check if the user has specified ystart & ystop
if ~isempty(varargin) && isa(varargin{1},'numeric') && length(varargin{1})==2
    ystart = varargin{1}(1);
    ystop  = varargin{1}(2);
    varargin(1) = [];
else
    ystart = min(y);
    ystop  = max(y);
end

% parse the remaining input properties specified by the user
for i = 1:length(varargin)
    if isa(varargin{i},'char')
        property_name = lower(varargin{i});
        if strcmp(property_name, 'labels')
            XYZlabels = varargin{i+1};
        elseif strcmp(property_name, 'colors')
            linecolors = varargin{i+1};
        elseif strcmp(property_name, 'ticklength')
            newTickLength = varargin{i+1};
        elseif strcmp(property_name, 'tickstyle')
            tickstyle = varargin{i+1};
        elseif strcmp(property_name, 'topmargin')
            axMup = varargin{i+1};
        elseif strcmp(property_name, 'bottommargin')
            axMdown = varargin{i+1};
        elseif strcmp(property_name, 'subplotgap')
            axMgap = varargin{i+1};
        elseif strcmp(property_name, 'deltaz')
            deltaz = varargin{i+1};
        elseif strcmp(property_name, 'subtractmean')
            subtractmean = true;
        end
    end
end

if ~exist('XYZlabels','var')
    XYZlabels = {'x','y','z'};
end

if ~exist('linecolors','var')
    linecolors = lines(Nslices);
end

if ~exist('newTickLength','var')
    newTickLength = 0.02;
end

if ~exist('tickstyle','var')
    tickstyle = '-';
end

if ~exist('axMup','var')
    axMup = 0.05; % upper margin of slice subplots
end

if ~exist('axMdown','var')
    axMdown = 0.1; % bottom margin of slice subplots
end

if ~exist('axMgap','var')
    axMgap = 0.005; % vertical gap between subplots
end

if exist('deltaz','var')
    axis_handles = gobjects(1,2);
    if ~exist('subtractmean','var')
        subtractmean = false;
    end
else
    axis_handles = gobjects(1,Nslices+1);
end

[~,ystart] = min(abs(y-ystart));
[~,ystop] = min(abs(y-ystop));
yslices = fliplr(round(linspace(ystart,ystop,Nslices)));

ax0L = 0.1; % left edge of image (normalized to figure dimensions)
ax0B = 0.25; % bottom edge of image (normalized to figure dimensions)
ax0W = 0.375; % image width (normalized to figure dimensions)
ax0H = 0.5; % image height (normalized to figure dimensions)
Hpixel = ax0H/length(y); % height of 1 pixel in image (normalized to figure dimensions)

figure;
axis_handles(1) = axes('position',[ax0L, ax0B, ax0W, ax0H]);
imagesc(x, y, z)
hold on
axis xy tight
set(gca,'fontsize',12,'linewidth',1)
xlabel(XYZlabels{1})
ylabel(XYZlabels{2})
set(gcf,'color','w')

hcb = colorbar('location','northoutside');
set(axis_handles(1),'position',[ax0L, ax0B, ax0W, ax0H]);
set(hcb,'fontsize',12,'linewidth',1)
xlabel(hcb, XYZlabels{3})

fpos = get(gcf,'position');
set(gcf,'position',[fpos(1:2)-fpos(3:4)*0.5, fpos(3:4)*1.5])

axLeft = 0.55; % axis left edge
axRight = 0.1; %axis right margin

axW = 1 -axRight -axLeft; % axis width
connection_endpoints = zeros(1,Nslices);

if ~exist('deltaz','var')
    axH = (1 -axMup -axMdown -(Nslices-1)*axMgap)/Nslices; % axis height
    axY1 = axMdown +(Nslices-1)*(axH +axMgap); % bottom edge of 1st axis
    for i = 1:Nslices
        axYi = axY1 -(axH +axMgap)*(i-1);
        connection_endpoints(i) = axYi +axH*0.5;   
        axis_handles(1+i) = axes('position', [axLeft, axYi, axW, axH]);
        plot(x, z(yslices(i),:), 'linewidth', 1, 'color', linecolors(i,:))
        
        set(gca,'fontsize',12,'linewidth',1,'yaxislocation','right')
        axis([min(x), max(x), min(min(z)), max(max(z))])
        box on
        if i < Nslices
            set(gca,'xticklabels','')
        else
            xlabel(XYZlabels{1})
            hyl = ylabel(XYZlabels{3});
            set(hyl,'units','normalized')
            hylpos = get(hyl,'position');
            hylDeltaY = (Nslices-1)*(axMgap+axH)/(2*axH);
            set(hyl,'position',[hylpos(1), hylpos(2)+hylDeltaY, hylpos(3)])
        end
        
        if mod(Nslices-i,2) == 1 && Nslices > 3
            set(gca,'yticklabels','')
        end
    end
    
else
    axH = 1 -axMdown -axMup;
    axis_handles(2) = axes('position', [axLeft, axMdown, axW, axH]);
    hold on
    for i = 1:Nslices
        z0 = deltaz*(Nslices-i);
        zi = z(yslices(i),:);
        if subtractmean
            zi = zi -mean(zi);
        end
        plot(x, z0+zi, 'linewidth', 1, 'color', linecolors(i,:))
    end
    set(gca,'fontsize',12,'linewidth',1,'yaxislocation','right')
    xlabel(XYZlabels{1})
    ylabel(XYZlabels{3})
    box on
    
    if subtractmean
        ax = [min(x), max(x), min(min(z')-mean(z')), ...
                deltaz*(Nslices-1)+max(max(z')-mean(z'))];
    else
        ax = [min(x), max(x), min(min(z)), deltaz*(Nslices-1)+max(max(z))];
    end
    
    axis(ax);
    for i = 1:Nslices
        z0 = deltaz*(Nslices-i);
        connection_endpoints(i) = axMdown +axH*(z0-ax(3))/(ax(4)-ax(3));
    end
end

% create the colored ticks & connecting lines
for i = 1:Nslices
    x1norm = ax0L +ax0W;
    y1norm = ax0B +Hpixel/2 +(ax0H-Hpixel)*(y(yslices(i))-min(y))/(max(y)-min(y));
    x2norm = axLeft;
    y2norm = connection_endpoints(i);
    annotation('line', [x1norm, x2norm], [y1norm, y2norm], 'linewidth', 1, 'color', linecolors(i,:))
    if newTickLength > 0
        annotation('line', [x1norm, x1norm-newTickLength*ax0W], ...
            y1norm*[1, 1], 'linewidth', 1, 'color', linecolors(i,:), ...
            'linestyle', tickstyle);
        if exist('deltaz','var')
            annotation('line', [x2norm, x2norm+0.02*ax0W], ...
                y2norm*[1, 1], 'linewidth', 1, 'color', linecolors(i,:), ...
                'linestyle', tickstyle);
        end
    end
end

end