
% ----------------------------------------------------------------------- %
% rho d multimodel class
% ----------------------------------------------------------------------- %
% This class contains a combined TBD model.
%
%
% Last update by Pablo: 07/2018 (maybe not accurate).
%
% P. I. Morales Guzman
% TUM / MPP
%
% ----------------------------------------------------------------------- %

classdef RHOD < handle
    
    properties (Constant = true, Hidden = true)
        
    end
    
    properties (Access=protected)
        
    end
    
    properties (Dependent=true,Access=public)
        
    end
    
    properties (Access=public)
        TBDCell;
        
        TBDDS;
        TBDIS;
        
        qUfrac;
        TimeSec;
        
        FPD_Segmentation;
        
        Q;
        Q_i = 18575;
        BKG_RateSec_i;
        
        nPixels;
        nqU;
        
        nModels;
        
        DTNormGS_i = 0;
        DTNormES_i = 0;
        
    end
    
    methods
        function obj = RHOD(varargin)
            %-------------------Parser Start------------------------------%
            
            p = inputParser;
            % Main inputs (necessary)
            p.addParameter('TBDCell',{},@(x) iscell(x)); %Study Object
            
            p.parse(varargin{:});
            
            obj.TBDCell              = p.Results.TBDCell;
            
            %-------------------Parser End--------------------------------%
            
            obj.nModels = length(obj.TBDCell);
            obj.nPixels = 4;
            obj.nqU = obj.TBDCell{1}.nqU;
            for mm = 1:obj.nModels
                obj.BKG_RateSec_i(1,mm) = obj.TBDCell{mm}.BKG_RateSec_i;
                obj.qUfrac(:,mm) = obj.TBDCell{mm}.qUfrac;
                obj.TimeSec(1,mm) = obj.TBDCell{mm}.TimeSec;
                
            end

            
        end % constructor     
    end % methods
    
    methods
        
        function ComputeTBDDS(obj,varargin)
            for mm = 1:obj.nModels
                % Inputs
                p = inputParser;
                p.addParameter('mSq_bias',0,@(x)isfloat(x));
                p.addParameter('E0_bias',0,@(x)isfloat(x));
                p.addParameter('N_bias',0,@(x)isfloat(x));
                %                 p.addParameter('DE_bias',0,@(x)isfloat(x));
                p.addParameter('B_bias',0,@(x)isfloat(x));
                %                 p.addParameter('TTGS_bias',0,@(x)isfloat(x));
                %                 p.addParameter('TTES_bias',0,@(x)isfloat(x));
                p.addParameter('DTGS_bias',0,@(x)isfloat(x));
                p.addParameter('DTES_bias',0,@(x)isfloat(x));
                %                 p.addParameter('mnu4Sq_Bias',0,@(x)isfloat(x));
                %                 p.addParameter('sin2T4_Bias',0,@(x)isfloat(x));
                %
                p.parse(varargin{:});
                mSq_bias = p.Results.mSq_bias;
                E0_bias  = p.Results.E0_bias;
                N_bias   = p.Results.N_bias;
                %                 obj.DE_Bias     = p.Results.DE_bias;
                B_bias   = p.Results.B_bias;
                %                 obj.TTGS_bias   = p.Results.TTGS_bias;
                %                 obj.TTES_bias   = p.Results.TTES_bias;
                %                 obj.DTGS_bias   = p.Results.DTGS_bias;
                %                 obj.DTES_bias   = p.Results.DTES_bias;
                %                 obj.mnu4Sq_Bias = p.Results.mnu4Sq_Bias;
                %                 obj.sin2T4_Bias = p.Results.sin2T4_Bias;
                
                obj.TBDCell{mm}.ComputeTBDDS(...
                    'mSq_bias',mSq_bias,'E0_bias',E0_bias,...
                    'N_bias',N_bias(mm),'B_bias',B_bias(mm));
            end
        end
        
        function ComputeTBDIS(obj)
            for mm = 1:obj.nModels
                obj.TBDCell{mm}.ComputeTBDIS();
                obj.TBDIS(:,mm) = obj.TBDCell{mm}.TBDIS;
            end
            obj.Q = obj.TBDCell{1}.Q;
            
        end
        
        
    end
end % class
