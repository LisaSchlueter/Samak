function [qU, hvdatamap] = read_hvdata(varargin)
%
% Read Krypton HV Data 
%
% Input   :
%  - Display (ON/OFF) 
%  - MK35    (Value)    - HV Calibration Constant
%
% Output  :
%  - qU                 - 1-D Array
%                         * Retarding potential in V
%  - hvdatamap          - 3-D Array
%                         * counts       in cps
%                         * count errors in cps
%                         * retarding potential accounting for 
%                           analysis plane + misaling correction
%                           (pixel-wise - from C. Weinheimer)
%
% Th. Lasserre - CEA Saclay
% Last Update: January 2018
%

% Parser
p = inputParser;
p.addParameter('Display','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('MK35',1972.4531);
p.parse(varargin{:});
Display    =    p.Results.Display;
MK35       =    p.Results.MK35;

% File Name
folder = 'hvdata';
prefix = 'KryptonRun-33149_pixel-';

for i=0:1:39
    Pixel = num2str(i);
    hvfile = sprintf('./%s/%s%s.txt',folder,prefix,Pixel);
    switch Display
        case 'ON'
            fprintf(2,'--- Opening %s ...\n',hvfile);
            figure(1)
    end
    hvdata = sortrows(importdata(hvfile),1);
    
    % read qU
    if (i==0)
        % Init
        hvdatamap = zeros(40,2,numel(hvdata(:,1)));
        qU        = flip(-(hvdata(:,1)) * MK35);
        switch Display
            case 'ON'
                fprintf(2,'--- qU begin\n');
                disp(qU');
                fprintf(2,'---qU end \n');
        end
    end
    
    % read pixel
    hvdatamap(i+1,1,:)    = flip(squeeze((hvdata(:,2))));
    hvdatamap(i+1,2,:)    = flip(squeeze((hvdata(:,3))));
    hvdatamap(i+1,3,:)    = flip(-(hvdata(:,1)) * MK35);
    
    switch Display
        case 'ON'
            fprintf(2,'--- Pixel %s begin\n',Pixel);
            hold on
            c=squeeze(hvdatamap(i+1,1,:));
            d=squeeze(hvdatamap(i+1,3,:));
            disp((c)');
            plot(d,c);
            fprintf(2,'--- Pixel %s end \n',Pixel);
    end
end

% Close figure
switch Display
    case 'ON'
        xlabel('qU');
        ylabel('counts per second');
        grid on
        title('L3-32 Line - Run 33149 - HV Calibration Paper - 3 Inner Rings');
        hold off
end

end

