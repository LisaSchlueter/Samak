
% ----------------------------------------------------------------------- %
% Class describing the KATRIN Experiment and Time Distribution (TD)
% ----------------------------------------------------------------------- %
% Building, Reading and Displaying TDs
%
% 3 Modes: Sim, DataKr, DataTBD, Read
% Input needed for Data: TD
% Input needed for Sim: qUmin, qUmax, qUStepSize, TDProfile
% Input needed for Read: TD-Name from simulation/katrinsetup/TD_DataBank
% If you want to simulate data: Use 'Mode' Data for the same TD
%
% Lisa Schlueter
% MPP/TUM
% March-2018
%
% Pablo I. Morales G.
% MPP/TUM
% 07/2018
%
%
% Th. Lasserre
% CEA/Saclay - TUM - IAS - MPP
%
% ----------------------------------------------------------------------- %

classdef KATRIN < handle %!dont change superclass without modifying parsing rules!
    
    properties (Constant = true, Hidden = true)
        % Physics Constants
        Year2Sec = 31557600;
    end
    
    properties (Access=protected)
        
    end
    
    properties (Dependent=true,Access=public)
        
    end
    
    properties (Access=public)
        %General Settings
        TDMode;           % DataKr, DataTBD, Sim, Read
        TDProfile;      % flat,....
        TD;             % Kr: Name of TD  &  TBD: Runnumber, for Sim: no meaning
        
        % Kinematics
        qU;             % retarding potential in V
        nqU;            % Number qU-Bins
        qUmin;          %
        qUmax;          %
        qUStepSize;     % qU-Bin-Width
        qUfrac;         % Time fraction spent at each qU
        
        % Generic
        Normalization;
        
        % Data Taking
        TimeSec;
        CPS;            % Output of IS in counts per second
        % Output mode
        quiet;
    end
    
    methods
        function obj = KATRIN(varargin)
            p = inputParser;
            p.addParameter('TDMode', 'Read', @(x)ismember(x,{'DataKr', 'DataTBD', 'Sim', 'Read','DScomparison','RFcomparison'}));
            p.addParameter('TD','Flat30'); %For DataKr: TD-Name; For DataTBD: Runnumber ; For Sim: meaningless
            p.addParameter('Normalization','NOMINAL',@(x)ismember(x,{'NOMINAL','COUNTS'}));
            p.addParameter('TimeSec',94672800,@(x)isfloat(x)); %3y
            p.addParameter('CPS','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('qU','',@(x)isfloat(x) && all(x>0,'all')); % if empty, read from data files
            p.addParameter('qUfrac','',@(x)isfloat(x));
            p.addParameter('qUmin',18545,@(x)isfloat(x) && x>0);
            p.addParameter('qUmax',18580,@(x)isfloat(x) && x>0);
            p.addParameter('qUStepSize',1,@(x)isfloat(x) && x>0);
            p.addParameter('TDProfile', 'flat', @(x)ismember(x,{'flat'}));
            p.addParameter('quiet','ON',@(x)ismember(x,{'ON','OFF'}));  % Output mode
            
            p.parse(varargin{:});
            
            obj.TDMode          = p.Results.TDMode;
            obj.Normalization   = p.Results.Normalization;
            obj.TimeSec         = p.Results.TimeSec;
            obj.CPS             = p.Results.CPS;
            obj.TD              = p.Results.TD;
            obj.qUmin           = p.Results.qUmin;
            obj.qUmax           = p.Results.qUmax;
            obj.qUStepSize      = p.Results.qUStepSize;
            obj.TDProfile       = p.Results.TDProfile;
            obj.quiet           = p.Results.quiet;
            obj.qU              = p.Results.qU;
            obj.qUfrac          = p.Results.qUfrac;
            if strcmp(obj.quiet, 'OFF'); fprintf(2,'Finished KATRIN.m constructor.\n'); end
            
        end % constructor
        
    end % methods
    
    methods
        
        function InitializeTD(obj,varargin)
            if ~isempty(obj.qU)
                %if MTD is given (qU and qUfrac) as input: don't read from file
                obj.qUmin      = min(min(obj.qU));
                obj.qUmax      = max(max(obj.qU));
                obj.nqU        = size(obj.qU,1);
                obj.qUStepSize = (obj.qU(end)-obj.qU(1)+1)/obj.nqU; 
                return;
            end
            
            switch obj.TDMode
                case 'DataKr'  %Krypton Data
                    % Change default TD name
                    if strcmp(obj.TD,'Flat30')
                        obj.TD = 'HVL332';
                        fprintf('Changing TD from Flat30 to HVL332 for Kr studies \n');
                    end
                    %   Read Data
                    if strcmp(obj.TD,'HVL332')
                        try   [qu, ~]= read_hvdata('Display','OFF','MK35',1972.4531);
                        catch; fprintf('ERROR: INVALID TD! \n');  end
                    else
                        try [qu, ~] = readKrData('TD',obj.TD);
                        catch; fprintf('ERROR: INVALID TD! \n'); end
                    end
                    
                    %Fill Parameters
                    obj.qU = qu;
                    td = ones(length(obj.qU),1); %all flat
                    obj.qUfrac = td/sum(td);
                    obj.qUmin = min(obj.qU); obj.qUmax=max(obj.qU);
                    obj.qUStepSize = (obj.qU(end)-obj.qU(1))/length(obj.qU); %mean Step Size
                    obj.nqU = length(obj.qU);
                    
                case 'DataTBD'    % Tritium Data
                    % everything read from file            
                case 'Sim'
                    obj.TD = 'Simulation';
                    qUSteps = ceil((obj.qUmax-obj.qUmin)/obj.qUStepSize)+1; %number of qU-Bins
                    obj.qU = linspace(obj.qUmin, obj.qUmax, qUSteps)';
                    obj.qUStepSize = (obj.qUmax-obj.qUmin + 1)/qUSteps; %new step size (<= init StepSize)
                    obj.nqU = length(obj.qU);
                    switch obj.TDProfile
                        case 'flat'
                            td = ones(length(obj.qU),1);
                            obj.qUfrac = td/sum(td);
                        case 'first_tritium'
                            % write first tritium TD
                    end
                case 'Read'
                    try  ReadTD(obj)
                    catch
                        fprintf(2,'KATRIN:InitializeTD: INVALID TD! \n');
                        return
                    end
            end
        end
        
        function ReadTD(obj)
            if      contains(obj.TD,'MTDcreator') || ...
                    contains(obj.TD,'Flat') || ...
                    contains(obj.TD,'DR30') || ...
                    contains(obj.TD,'Sensitivity')  || ...
                    contains(obj.TD,'Optimize') || ...
                    contains(obj.TD,'KITNuScan') || ...
                    contains(obj.TD,'MasterRF') || ...
                    contains(obj.TD,'KNM1_25p') 
                file = [obj.TD,'.mat'];
            else
                if contains(obj.TD,'Fake')
                    file = sprintf('%s.mat', obj.TD);
                else
                    file = sprintf('Run%s.mat', obj.TD);
                end
            end
            filename = [getenv('SamakPath'),'/simulation/katrinsetup/TD_DataBank/',file];
            myTD = load(filename);
            obj.qU = myTD.qU;
            obj.qUfrac= myTD.qUfrac;
            obj.qUmin = min(min(obj.qU)); obj.qUmax = max(max(obj.qU));
            obj.qUStepSize = (obj.qU(end)-obj.qU(1))/length(obj.qU); %mean Step Size
            obj.nqU = size(obj.qU,1);
        end
        
        function PlotTD(obj,varargin)
            p = inputParser;
            p.addParameter('ring',1,@(x)isfloat(x));
            p.addParameter('xLim','',@(x)isfloat(x));
            p.addParameter('yLim','',@(x)isfloat(x));
            p.addParameter('Color','Amethyst',@(x)ischar(x));
            p.parse(varargin{:});
            ring      = p.Results.ring;
            Color = p.Results.Color;
            xLim = p.Results.xLim;
            yLim = p.Results.yLim;
            f1 = figure('Name','MTD','Renderer','opengl');
            set(f1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
            bar(obj.qU(:,ring)-18574, obj.qUfrac(:,ring),'FaceColor',rgb(Color));
           PrettyFigureFormat('FontSize',22);
            %title(sprintf('MTD: %s',strrep(obj.TD,'_',' ')),'FontWeight','normal');
            leg = legend(sprintf('MTD: %s',strrep(obj.TD,'_',' ')),'Location','northeast');
            PrettyLegendFormat(leg);
           % leg.FontSize = get(gca,'FontSize');
            xlabel(sprintf('Retarding potential qU - %.0f (eV)',18574));
            ylabel('Time fraction');

         
         
            if ~isempty(xLim)
                xlim(xLim);
            end
            
             if ~isempty(yLim)
                ylim(yLim);
             else
                    ylim([0.000 max(max(obj.qUfrac))*1.1]);
             end
            % save 
               ptlfile = sprintf('MTD-%s',obj.TD);
          
            if ~exist('./plots/MTD/','dir')
                mkdir ./plots/MTD
            end
            %publish_figurePDF(f1,['plots/MTD/pdf/',ptlfile,'.pdf']);
            print(f1,['./plots/MTD/',ptlfile,'.png'],'-dpng','-r250');
            fprintf('save plot to %s \n',['./plots/MTD/',ptlfile,'.png']);
        end
        
    end
end % class
