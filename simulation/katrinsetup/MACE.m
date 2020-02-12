% ----------------------------------------------------------------------- %
% Class describing the KATRIN Magnetic Adiabatic Cooling Filter (MACE)
% ----------------------------------------------------------------------- %
% 
% - MACE parameters
% - Transmission Function
% - KATRIN Response Function
% - Covariance Matrix for WGTS+MACE
% 
% 
% Th. Lasserre, 2017
% CEA - Saclay
% TUM - IAS
%
% ----------------------------------------------------------------------- %

%classdef MACE < WGTS & TBD
classdef MACE < handle
   
    properties (Constant = true, Hidden = true)
        
        % Physics Constants

    end
    
    properties (Access=protected)
        
        % Mac-E Filter
        Bm       = 6;              % Pinch B field, Tesla
        Bs       = 3.6;            % Source B field, Tesla
        Ba       = 0.0003;         % Analysis Plane B field, Tesla
        MacEres  = 0.93;           % Mac-E Energy resolution, eV
        
    end
    
    properties (Dependent=true,Access=public)
    end
    
    properties (Access=public)
        yyy;
    end
    
    methods
        function obj = MACE(varargin)
            fprintf(2,'Processing MACE Constructor ...\n');
            p = inputParser;
            
            % Parameters
            p.addParameter('yyy',0,@(x)isfloat(x) && x>0);

            p.parse(varargin{:});
            
            % Parameters
            obj.yyy  = p.Results.yyy;
            
        end % constructor
        
    end % methods
    
    methods
        
    end
end % class
