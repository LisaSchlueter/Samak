function [plotHandle, cbHandle] = FPDViewer(data,varargin)
%%
%
% P. I. Morales, December 2017
%

% datasize = 148 or 13

p = inputParser;

p.addParameter('ReDrawSkeleton','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('outliers',-1,@(x)isfloat(x) && min(x)>=0);

p.parse(varargin{:});

ReDrawSkeleton  = p.Results.ReDrawSkeleton;
outliers        = p.Results.outliers;

if strcmp(ReDrawSkeleton,'OFF')
    fig = openfig('FPDFrame.fig');
end
fig = get(groot,'CurrentFigure');

%% Parameters of the detector
% Start and end of radius [cm]
rStart = [0 0.7398 1.4796 1.9573 2.3394 2.6674 2.9592 3.2247 ...
    3.4699 3.699 3.9146 4.119 4.3137]';
rEnd = [0.7398 1.4796 1.9573 2.3394 2.6674 2.9592 3.2247 ...
    3.4699 3.699 3.9146 4.119 4.3137 4.5]';
% number of segments per ring
segs = [4,12*ones(1,12)]';
% starting angle [rad]
thetaStart = ([0 -15 0 -15 0 -15 0 -15 0 -15 0 -15 0].*pi/180)';

%% Plot detector outline
if strcmp(ReDrawSkeleton,'ON')
    hold on
    % set the axes in the correct position
    % (important for the numbering in the detectors)
    set(gca,'XLim',[-max(rEnd) max(rEnd)], ...
        'YLim',[-max(rEnd) max(rEnd)],...
        'Position',[0.11452380952381 0.11 0.682738095238095 0.815]);
    pix = 0; %number in each detector
    % build the skeleton segment by segment from the center outwards
    for ii = 1:13
        for jj = 0:segs(ii)-1
            wS = 2*pi/segs(ii);
            th1 = jj*wS + thetaStart(ii);
            th2 = (jj+1)*wS + thetaStart(ii);
            % function the get the values of r and theta, according to the
            % begining and ending radius and angles in each segment
            [theta,r] = getWedgeBorder(th1,th2,rStart(ii),rEnd(ii));
            p1 = polar(theta,r);
            % xx and yy are the position of the numbers in each segment
            %         xx = 0.13+0.775/2 + ((rStart(ii)+rEnd(ii))/2*cos((th1+th2)/2))/(9/0.775);
            xx = 0.114523809523810+0.682738095238095/2 + ((rStart(ii)+rEnd(ii))/2*cos((th1+th2)/2))/(9/0.682738095238095);
            yy = 0.11+0.815/2 + ((rStart(ii)+rEnd(ii))/2*sin((th1+th2)/2))/(9/0.815);
            annotation(gcf,'textbox',[xx yy 0.0 0.0],'Color','k',...
                'String',num2str(pix),'FontSize',11,'Margin',0,...
                'HorizontalAlignment','center','VerticalAlignment','middle');
            pix = pix + 1;
        end
    end
 
    % Set all lines to color black
    tmp = findall(gca,'Type','line');
    detLines = tmp(1:148);
    %set(detLines,'Color','k');
    set(detLines,'Color',rgb('SteelBlue'),'LineWidth',2)
    hold off
end

%% Pixel or rings values
% Each segment or ring is filled according to the values given

% Check if data has column form
if iscolumn(data)
    data = data';
end

% Give outliers a NaN value so they appear white
if sum(outliers) >= 0
   % data(ismember(1:length(data),outliers+1)) = NaN;
%     outliers_count = 1;
%     for ii = 1:length(data)
%         if 
%         if outliers(outliers_count) == ii
%             data(ii,1) = NaN;
%             outliers_count = outliers_count + 1;
%         end
%         
%     end
end

% Fill the pixels in the FPD
hold on
if length(data) == 148 % for pixels
    
    pixels1 = 0; pixels2 = 4;
    
     for ii = 1:13
        fillBullseye(data(1+pixels1:pixels2),rStart(ii),rEnd(ii),thetaStart(ii),rgb('SlateGray'));
        pixels1 = pixels2;
        pixels2 = pixels2 + 12;
     end
     
     pixels1 = 0; pixels2 = 4;
   
    for ii = 1:13
        fillBullseye(data(1+pixels1:pixels2),rStart(ii),rEnd(ii),thetaStart(ii));
        pixels1 = pixels2;
        pixels2 = pixels2 + 12;
    end
elseif length(data) == 13 % for rings
    for ii = 1:13
        fillBullseye(data(ii),rStart(ii),rEnd(ii),thetaStart(ii));
    end
end

% Bring the skeleton of the detector to the top
if strcmp(ReDrawSkeleton,'OFF')
    tmp = findall(gca,'Type','line');
    detLines = tmp(1:148);
end
uistack(detLines,'top');
% Colorbar properties adjusted so that the numbers stay in place
cbHandle = colorbar;
cbHandle.Color = 'k'; % Change color according to convenience
cbHandle.Location = 'manual';
cbHandle.Position = [0.8333 0.1413 0.0476 0.7524];
plotHandle = gcf;
% new: style
PrettyFigureFormat;
box off; % no upper box around FPD
set(gca,'Visible','off'); % no axes
plotHandle.Units = 'normalized';
plotHandle.Position = [0.1,0.1,0.42,0.55]; % make it rounder

hold off

end

function [theta_wedge,rho_wedge] = getWedgeBorder(theta1,theta2,rho1,rho2)
% make the points for the angle and the radius so that it looks smooth
dtheta = (theta2 - theta1)/(4*pi);
drho = (rho2 - rho1)/12;
arc = theta1:dtheta:theta2;
lin = rho1:drho:rho2;

% create the values used for the polar plot.
theta_wedge = [arc theta2*ones(1,13)...
    arc(end:-1:1) theta1*ones(1,13)];
rho_wedge = [ones(1,length(arc))*rho1 lin ...
    ones(1,length(arc))*rho2 lin(end:-1:1)];

end

function surfaceObject = fillBullseye(cdata,r1,r2,theta1,color)
segsNum = length(cdata);

switch segsNum
    case 1
        upsmp = 72;
    case 4
        upsmp = 18;
    case 12
        upsmp = 6;
end

%Upsample cdata
cdata = reshape(repmat(cdata,upsmp,1),1,segsNum*upsmp);
theta = repmat(theta1:5*pi/180:theta1+2*pi,2,1);
lin = [r2 r1];

% Create meshgrid for surface plot.
X = repmat(lin',1,segsNum*upsmp+1).*cos(theta);
Y = repmat(lin',1,segsNum*upsmp+1).*sin(theta);

surfaceObject = surf(gca,X,Y,zeros(size(X)),cdata);
set(surfaceObject,'EdgeColor','none');

if exist('color','var')
  set(surfaceObject,'FaceColor',color);  
end

end
