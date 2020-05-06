classdef SterileAnalysis < handle
    % class to collect and organise sterile analysis routines
    %%
    properties (Access=public)
        RunAnaObj; % Mother Object: defines RunList, Column Density, Stacking Cuts,....
    end
    
    methods % constructor
        function obj = SterileAnalysis(varargin)
            p = inputParser;
            p.addParameter('RunAnaObj','', @(x) isa(x,'RunAnalysis') || isa(x,'MultiRunAnalysis'));
            p.parse(varargin{:});
            obj.RunAnaObj = p.Results.RunAnaObj;
            
            if isempty(obj.RunAnaObj)
                fprintf('RunAnaObj has to be specified! \n');
                return
            end
            GetSamakPath; %sets current samak path as enviromental variable
        end  
    end
    
    methods 
        function GridSearch(varargin) %
            p = inputParser;
            p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            RecomputeFlag = p.Results.RecomputeFlag;
        end        
    end
    
    
end