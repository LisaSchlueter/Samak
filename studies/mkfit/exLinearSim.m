%% Examplary simulation model for calculating a linear graph. 
%
% Note: For mimimizations based on using '../fitting/fminuit.m and '../minimize.m' the simulation 
%       model must provide a method (here: 'y=compute(obj,varargin)') that accepts parameters in
%       pairs of ''name',value' and returns an array which reflects the simulated theory.
%
% 2017 Nov - Marc Korzeczek

classdef exLinearSim
    properties
       x; 
    end
    
    methods    
        function obj=exLinearSim()           
            %empty initializer
             obj.x= .1:.1:10;
        end
    
        function y=compute(obj,varargin)
            %return a linear curve for input x-values
            %and two parameters slope & offset
            p=inputParser();
            p.addParameter('slope',1);
            p.addParameter('offset',0);
            p.addParameter('poisson',false);
            p.parse(varargin{:});
            
            slope=p.Results.slope;
            offset=p.Results.offset;
            poisson=p.Results.poisson;
            
            y=slope* obj.x +offset;
            if poisson ; y=y+sqrt(abs(y)).*randn(1,numel(y)); end
        end
    end
end
