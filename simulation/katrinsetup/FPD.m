% ----------------------------------------------------------------------- %
% Class describing the KATRIN Focal Plane Detector(FPD)
% ----------------------------------------------------------------------- %
%
% - FPD parameters
% - Integral Electron Spectrum
%   - FPD Detector-wise
%   - FPD Ring-wise
%   - FPD Pixe-wise
% - Covariance Matrix for Detector Systematics
%
%
% Th. Lasserre
% CEA/Saclay - TUM - IAS - MPP
%
% Last Updated: Jan. 2018
%
% ----------------------------------------------------------------------- %

classdef FPD < KATRIN & handle %!dont change superclass without modifying parsing rules!
    
    properties (Constant = true, Hidden = true)
        
    end
    
    properties (Access=protected)
        
    end
    
    properties (Dependent=true,Access=public)
    end
    
    properties (Access=public)
        
        % FPD: segmentation 
        FPD_Segmentation;        % Segmentation: OFF/RING/SINGLEPIXEL/MULTIPIXEL
        FPD_PixList;             % list of pixels considered (for normalization of activity, response function,...)
        FPD_RingList;            % list of rings considered 
        FPD_RingPixList;         % (not an input) FPD Ring <--> Pixels (cell)
        nPixels;                 % (not an input) number of pixels analyzed in multipixel analysis 
        nRings;                  % (not an input) number of rings analyzed in multipixel analysis 
        FPD_RingMerge;           % How the rings are merged
        
        % other FPD attributes
        FPD_MeanEff          ;        % Detector Efficiency
        FPD_Coverage         ;        % Coverage of Detector: depends on ROI cut
        FPD_ROIlow           ;        % Region of Interest (ROI)- lower cut [keV] -> affects efficiency
        FPD_ROIup=32         ;        % Region of Interest (ROI)- upper cut [keV]
        FPD_ROIEff           ;        % Correction for ROI qU-dependent efficiency
        FPD_PileUpEff        ;        % Correction for PileUp pixel rate-dependent efficiency
        FPD_Eff_pix          ;        % Detector Efficiency per Pixel, vector
         
        % Backgrounds
        BKG_Flag            ;         % ON/OFF
        BKG_Type            ;         % Flat/Slope
        BKG_RateAllFPDSec   ;         % Rate per second per FPD
        BKG_RateRingSec     ;         % Rate per Second per ring
        BKG_RatePixelSec    ;         % Rate per second per pixel
        BKG_RateSec;  BKG_RateSec_i;  % Used for computations
        BKG_RateSecallPixels; BKG_RateSecallPixels_i; % Pixel-wise rate per second
        BKG_RateSecErrallPixels;      % Pixel-wise rate per second error
        BKG_Slope;  BKG_Slope_i=0;
    end
    
    methods
        function obj = FPD(varargin)
            p = inputParser;
            % FPD Parameters
            p.addParameter('FPD_Segmentation','OFF',@(x)ismember(x,{'OFF','RING','MULTIPIXEL','SINGLEPIXEL'}));
            p.addParameter('FPD_PixList',1:148,@(x)isfloat(x) && all(x)>0);
            p.addParameter('FPD_RingList',1:13,@(x)isfloat(x) && all(x)>0);
            p.addParameter('FPD_MeanEff',0.95,@(x)isfloat(x) && x>0);
            p.addParameter('FPD_Coverage',1,@(x)isfloat(x) && x>0);
            p.addParameter('FPD_Eff_pix',0,@(x)isfloat(x) && sum(x)>0 && sum(x)<=148);
            p.addParameter('FPD_ROIEff','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('FPD_PileUpEff','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.addParameter('FPD_ROIlow','',@(x)isfloat(x) && x>0); %if empty: do nothing, if not: compute coverage 
            p.addParameter('FPD_RingMerge','Default',@(x)ismember(x,{'Default','None','Full','Half','Azi','AziHalf'}));
    
            % Background Parameters
            p.addParameter('BKG_Flag','ON',@(x)ismember(x,{'ON','OFF','XmasData'}));
            p.addParameter('BKG_Type','FLAT',@(x)ismember(x,{'FLAT','SLOPE'}));
            p.addParameter('BKG_RateAllFPDSec',300e-3,@(x)all(isfloat(x)) && x>=0);
            p.addParameter('BKG_RateRingSec',10e-3,@(x)all(isfloat(x)));
            p.addParameter('BKG_RatePixelSec',1e-3,@(x)all(isfloat(x)));
            
            %|<- Inherited KATRIN.m ...
            %p.addParameter('Normalization','NOMINAL',@(x)ismember(x,{'NOMINAL','COUNTS'}));
            %p.addParameter('TimeYear',3,@(x)isfloat(x) && x>0);
            %p.addParameter('TD','DR30',@(x)ismember(x,{'DR20','DR30','DR40','DR50','Kle15','Kle15Ext','Kle15Dyn','Flat20','Flat30','Flat50','Flat60','Flat100'}));
            %p.addParameter('quiet','OFF',@(x)ismember(x,{'ON','OFF'}));
            %... Inherited KATRIN.m ->|
            
            % Parse unmatched parameters to KATRIN.m
            p.KeepUnmatched=1;
            p.parse(varargin{:});
            if( isempty(fieldnames(p.Unmatched)) ); Unmatched={}; else
                Unmatched = reshape(...
                    [fieldnames(p.Unmatched),struct2cell(p.Unmatched)]',...
                    [1,length(fieldnames(p.Unmatched))*2]);
            end
            obj=obj@KATRIN(Unmatched{:}); %Parse to superclass KATRIN.m
            
            % FPD Parameters
            obj.FPD_Segmentation    = p.Results.FPD_Segmentation;
            obj.FPD_PixList         = p.Results.FPD_PixList;
            obj.FPD_RingList        = p.Results.FPD_RingList;
            obj.FPD_MeanEff         = p.Results.FPD_MeanEff;
            obj.FPD_Coverage        = p.Results.FPD_Coverage;
            obj.FPD_ROIlow          = p.Results.FPD_ROIlow;
            obj.FPD_Eff_pix         = p.Results.FPD_Eff_pix;
            obj.FPD_ROIEff          = p.Results.FPD_ROIEff;
            obj.FPD_PileUpEff       = p.Results.FPD_PileUpEff;
            obj.FPD_RingMerge       = p.Results.FPD_RingMerge;
            
            % Background Parameters
            obj.BKG_Flag            = p.Results.BKG_Flag;
            obj.BKG_Type            = p.Results.BKG_Type;
            obj.BKG_RateAllFPDSec   = p.Results.BKG_RateAllFPDSec;
            obj.BKG_RateRingSec     = p.Results.BKG_RateRingSec;
            obj.BKG_RatePixelSec    = p.Results.BKG_RatePixelSec;
            
            % Synchronize/overlap selected pixels with selected rings
            switch obj.FPD_RingMerge
                case 'Default' % 10 rings
                    [obj.FPD_PixList,obj.FPD_RingPixList] = Ring2PixelDefCombi(obj.FPD_RingList,obj.FPD_PixList);
                case 'None' % 12 rings
                    [obj.FPD_PixList,obj.FPD_RingPixList] = Ring2Pixel(obj.FPD_RingList,obj.FPD_PixList);
                case 'Full' % 4 pseudo rings
                    [obj.FPD_PixList,obj.FPD_RingPixList] = Ring2PixelCombi(obj.FPD_RingList,obj.FPD_PixList);
                case 'Half'
                    [obj.FPD_PixList,obj.FPD_RingPixList] = Ring2PixelHalfCombi(obj.FPD_RingList,obj.FPD_PixList);
            end

            obj.nRings = numel(obj.FPD_RingPixList);
            switch obj.FPD_Segmentation
                case {'SINGLEPIXEL','MULTIPIXEL'}
            obj.nPixels = numel(obj.FPD_PixList);
                case 'OFF'
                    obj.nPixels = 1;
                case 'RING'
                    obj.nPixels = obj.nRings;
            end
            
            CreateqUmultipix(obj) % only used in simulation
            %SelectPixel(obj)

            if( strcmp(obj.quiet, 'OFF') ); fprintf(2,'Finished FPD.m constructor.\n'); end
            
            if ~isempty(obj.FPD_ROIlow)
                obj.FPD_Coverage =  GetFPDcoverage('FPD_ROIlow',obj.FPD_ROIlow);
            end
        end % constructor
        
    end % methods
    
    methods
        
        function CreateqUmultipix(obj)
            % If doing a multipixel simulation, extend qU to whole detector
            % (pixels are selected later)
            if strcmp(obj.TDMode,'Sim') && isvector(obj.qU)
                switch obj.FPD_Segmentation
                    % TODO: Add E and B corrections in the analysis plane
                    case 'MULTIPIXEL'
                        obj.qU = repmat(obj.qU,1,148);
                    case 'RING'
                        obj.qU = repmat(obj.qU,1,13);
                end
            end
        end
        
        function SelectPixel(obj)
            % Take qU values just for the pixels to be analyzed. 
            % In the case of 'OFF' (detector as one big pixel), the qU are
            % averaged.
            
            switch obj.FPD_Segmentation
                case 'OFF'
                    %fprintf('Calculating average of qU of all pixels. \n')
                    obj.qU = mean(obj.qU,2);
                    %fprintf('Changing pixels to pixel 1 (to make the rest of the calculations easier). \n')
                    obj.FPD_Pixel = 1;
                case 'RING'
                    % Work in progress
                    % For simplicity, it will only work with one ring for
                    % now.
                    % TODO: Make multiring analysis
                    obj.nRings = length(obj.FPD_RingList);
                    if obj.nRings > 1
                        error('You can only input one ring for the moment. Later updates may include the possibility to analyze more rings simultaneously.')
                    end
                    qUmpix = obj.qU;
                    obj.qU = zeros(obj.nqU,obj.nRings);
                    obj.InitFPDrings(); % this function gives obj.ring{1...13}
                    for ri = 1:obj.nRings
                        if size(qUmpix,2) > 1
                            obj.qU(:,ri) = mean(qUmpix(:,obj.ring{obj.FPD_RingList(ri)}),2);
                        else
                           obj.qU(:,ri) = qUmpix; 
                        end
                    end
                    obj.FPD_Pixel = 1; 
                case 'MULTIPIXEL'
                    obj.qU = obj.qU(:,obj.FPD_Pixel);
                case 'SINGLEPIXEL'
                    if ~isscalar(obj.FPD_Pixel)
                        warning('Changing pixels to pixel 1. If you want to analyze other pixel, indicate it in FPD_Pixel.')
                        obj.FPD_Pixel = 1;
                    end
            end
            obj.nPixels = length(obj.FPD_Pixel);
        end
        
        function InitializeDetectorEfficiency(obj)
            % If there is information on the efficiency per pixel, it is
            % taken into account here. It should be given as a vector in
            % 'obj.FPD_Eff_pix'
            
            switch obj.FPD_Segmentation
                case 'OFF'
                    if any(obj.FPD_Eff_pix ~= 0)
                        obj.FPD_MeanEff = mean(obj.FPD_Eff_pix);
                    end
                case 'RING'
                    % It takes the average efficiency of the pixels per
                    % ring.
                    if all(obj.FPD_Eff_pix == 0)
                        obj.FPD_Eff_pix = obj.FPD_MeanEff*ones(obj.nRings,1);
                    else
                        fprintf('Taking detector efficiency from vector.')
                        efficiencympix = obj.FPD_Eff_pix;
                        for ri = 1:obj.nRings
                            obj.FPD_Eff_pix(:,ri) = mean(efficiencympix(obj.ring{obj.FPD_RingList(ri)}));
                        end
                    end
                case {'SINGLEPIXEL','MULTIPIXEL'}
                    if all(obj.FPD_Eff_pix == 0)
                        obj.FPD_Eff_pix = obj.FPD_MeanEff*ones(obj.nPixels,1);
                    else
                        fprintf('Taking detector efficiency from vector.')
                        obj.FPD_Eff_pix = obj.FPD_Eff_pix(obj.FPD_Pixel);
                    end
            end
            
        end
            
        function InitializeBackground(obj,varargin)
            % Initializes background 
            % BCK chosen by the flag BKG_Flag: 
            %
            % 'XmasData' gives to the model
            % the background from the Chistmas runs 35023 to 35110; 
            %
            % 'ON' gives a flat background chosen in
            % 'BKG_RateSecallPixels_i' as cps;
            %
            % 'OFF' gives a flat background of 0.
            %
            % All backgrounds are then fluctuated according to a gaussian function 
            %
            
            % Parser
            p = inputParser;
            p.addParameter('Random','OFF',@(x)ismember(x,{'ON','OFF'}));
            p.parse(varargin{:});
            Random  = p.Results.Random;
            
            % 148 pixels hardcoded...
            obj.BKG_RateSecallPixels_i        = obj.BKG_RateAllFPDSec*ones(148,1)/148;
            
            switch Random
                case 'ON'
                    obj.BKG_RateSecallPixels_i  = obj.BKG_RateSecallPixels_i ...
                        + randn(1,148).* obj.BKG_RateSecErrallPixels;
            end
            
            switch obj.BKG_Flag
                case 'XmasData'
                    bck_struct = load([getenv('SamakPath'),'/inputs/BCK/BackgroundRate_35023-35110.mat']);
                    bck = bck_struct.bck;
                    bck(:,1) = bck(:,1) + 1;
                    bck(:,2) = normrnd(bck(:,2),bck(:,3));
                    if min(bck(:,2)) < 0
                        error('BACKGROUND NEGATIVE! STOP!');
                    end
                case 'ON'
                    bck(:,1) = (1:148)';
                    bck(:,2) = obj.BKG_RateSecallPixels_i;
                    bck(:,3) = sqrt(obj.BKG_RateSecallPixels_i);
                    if min(bck(:,2)) < 0
                        error('BACKGROUND NEGATIVE! STOP!');
                    end
                case 'OFF'
                    obj.BKG_RateSec_i = zeros(1,obj.nPixels);
                    return;
            end
            
            
            if strcmp(obj.BKG_Flag,'ON')
                switch obj.FPD_Segmentation
                    case 'OFF'
                        obj.BKG_RateSec_i = obj.BKG_RateAllFPDSec;
                    case 'RING'
                        obj.BKG_RateSec_i = obj.BKG_RateRingSec;
                    case {'MULTIPIXEL','SINGLEPIXEL'}
                        obj.BKG_RateSec_i = bck(obj.FPD_Pixel,2)';
                end
            elseif strcmp(obj.BKG_Flag,'XmasData')
                switch obj.FPD_Segmentation
                    case 'OFF'
                        obj.BKG_RateSec_i = sum(bck(:,2));
                    case 'RING'
                        obj.BKG_RateSec_i = cell2mat(cellfun(@(x) sum(bck(x)),obj.FPD_RingPixList,'UniformOutput',false)');
                    case {'MULTIPIXEL','SINGLEPIXEL'}
                        obj.BKG_RateSec_i = bck(obj.FPD_Pixel,2)';
                end
            end
        end

        function SetBkgRate(obj,value)
            switch obj.BKG_Flag
                case 'ON'
                    switch obj.FPD_Segmentation
                        case 'OFF'
                            obj.BKG_RateSec_i = value;
                        case 'RING'
                            obj.BKG_RateSec_i = value;
                        case 'PIXEL'
                            obj.BKG_RateSec_i = value;
                    end
                case 'OFF'
                    obj.BKG_RateSec_i = 0;
            end
        end
        
        
        function PlotBCK(obj,varargin)
            p = inputParser;
            p.addParameter('fign',999,@(x)isfloat(x) && x>0);
            p.addParameter('viewer','FPD',@(x)ismember(x,{'FPD','lin'}));
            
            p.parse(varargin{:});
            fign  = p.Results.fign;
            viewer= p.Results.viewer;
            
            switch obj.FPD_Segmentation
                case 'OFF'
                    figure(fign)
                    plot((1:148)',ones(1,148)*obj.BKG_RateSec_i);
                    xlabel('Pixel'); ylabel('background rate [cps]');
                case 'RING'
                    error('Not working yet. Try again in the next release.')
                case 'PIXEL'
                    if obj.FPD_Pixel == 0
                        switch viewer
                            case 'FPD'
                                FPDViewer(obj.BKG_RateSecallPixels_i)
                            case 'lin'
                                figure(fign)
                                plot(1:148,obj.BKG_RateSecallPixels_i);
                                xlabel('Pixel'); ylabel('background rate [cps]');
                        end
                    else
                        figure(fign)
                        plot(1:148,obj.BKG_RateSec_i);
                        xlabel('Pixel'); ylabel('background rate [cps]');
                    end
            end
            
        end %PlotBck
        
    end
    
    methods(Static)
                function ROIEffCorrectionFactor = qUEfficiencyCorrectionFactor(qu)
            %
            % qU dependent efficiency correction factor
            % Based on Sanshiro inputs [16-18.6] KV
            %
            % Input:
            % - qU (Volt)
            % Output
            % - Correction Factor (factorized with FPD efficiency)
            %
            %     Linear model Poly4:
            %      fitresult(x) = p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5
            %      Coefficients (with 95% confidence bounds):
            %        p1 =   -0.000114  (-0.0001155, -0.0001125)
            %        p2 =    0.007734  (0.007628, 0.007839)
            %        p3 =     -0.1967  (-0.1995, -0.194)
            %        p4 =       2.228  (2.196, 2.259)
            %        p5 =      -8.488  (-8.624, -8.352)
            %             %
            % Thierry Lasserre
            % Last Modfied: 27/07/2018
            %
            
            if min(qu) < 16.5
                return;
            end
            
            p1 =     -1.139949799271709e-04;
            p2 =      0.007733509590149;
            p3 =     -0.196747862772578;
            p4 =      2.227682510553492;
            p5 =     -8.488123971598990;
            
            x=qu/1000;
            ROIEffCorrectionFactor = ...
                p1*x.^4 + p2*x.^3 + p3*x.^2 + p4*x +p5;
        end
        
        function PUeffCorrectionFactor = PileUpEfficiencyCorrectionFactor(pixelrate)
            %
            % FPD efficiency correction factor due to Pile-Up
            % Depends on the rate per pixel
            % Based on Sanshiro inputs
            %
            % Input:
            % - pixel-rate (cps)
            % Output
            % - Correction Factor (factorized with FPD efficiency)
            %
            % Linear model Poly1:
            %      fitresult(x) = p1*pixelrate + p2
            %      Coefficients (with 95% confidence bounds):
            %        p1 =  -1.598e-06  (-1.598e-06, -1.597e-06)
            %        p2 =           1  (1, 1)
            %
            % Thierry Lasserre
            % Last Modfied: 27/07/2018
            %
            p1 =  -1.597653658180098e-06;
            p2 =  0.999999590625725;
            
            PUeffCorrectionFactor = ...
                p1.*pixelrate + p2;
        end
    end
end % class
