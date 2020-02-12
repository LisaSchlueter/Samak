%---------------------------------------------------------------------------------------------------
% Function to scale the Background (BKG) according to magnetic fields strength in
% the analyzing plane (Ba)
% Taken from Nikolaus Trost phD Thesis, computed for WGTS_B_T = 3.6 T
% Scaled b us to match any other WGTS_B_T field configuration
% Using: Bs * WGTS Section = Ba * Analysis plane section
% (reducing Bs means reducing Analysis plane section, and therefore any volumic background)
%
% ---- AnchorBkg6G Option:-----
% If AnchorBkg6G is empty:
% Background is scaled according to ROI cut, based on FT Background level
% If AnchorBkg6G is not empty:
% Background is scaled according AnchorBkg6G level
% ------------------------------------------------------------------------------------------
% Lisa Schlueter (September 2018)
% Inputs:
%  - WGTS_B_T
%  - MACE_Ba_T
%  - FPD_ROIlow or AnchorBkg6G (if both filled, AnchorBkg6G will be used)
% Output:
%  - Background in cps
%-------------------------------------------------------------------------------------------
function BKG = GetBackground(varargin)
% ------------------------------------------------------------------------------------------------
p = inputParser;

p.addParameter('MACE_Ba_T',6e-4,@(x)all(isfloat(x)) && all(x>=0));
p.addParameter('WGTS_B_T',3.6*0.7,@(x)all(isfloat(x)) && all(x>=0)); % B field reduction factor
%p.addParameter('MACE_Ba_T',[]);
%p.addParameter('WGTS_B_T',[]); % B field reduction factor
p.addParameter('AnchorBkg6G',''); % if empty, Anchor based und ROIlow is used 335=FT, [26keV, 32 keV]ROI, 238 Nikolaus PhD (70%bfields)
p.addParameter('FPD_ROIlow',14,@(x)isfloat(x));
p.addParameter('FPD_PixList',1:148,@(x)all(isfloat(x)) && all(x)>0);

p.parse(varargin{:});

MACE_Ba_T     = p.Results.MACE_Ba_T;
WGTS_B_T      = p.Results.WGTS_B_T;
AnchorBkg6G   = p.Results.AnchorBkg6G; % scales Background as such, that Bkg is @Ba=6G at this level
FPD_ROIlow    = p.Results.FPD_ROIlow;  % if AnchorBkg6G is empty -> uses ROI cut to get AnchorBkg6G based on FT level
FPD_PixList   = p.Results.FPD_PixList;  % Pixel List used for Unifrom Analysis

if numel(MACE_Ba_T)==1 && numel(WGTS_B_T)>=1            % In case Ba is scalar and BT is array
    BTsize = size(WGTS_B_T);
    MACE_Ba_T = repmat(MACE_Ba_T,BTsize(1),BTsize(2));
elseif numel(MACE_Ba_T)>=1 && numel(WGTS_B_T)==1        % In case BT is scalar and Ba is array
    BaSize = size(MACE_Ba_T);
    WGTS_B_T = repmat(WGTS_B_T,BaSize(1),BaSize(2));
elseif numel(WGTS_B_T)~=numel(MACE_Ba_T)                % In case both array with different size
    fprintf(2,'Error: Please insert Ba and B_T arrays of the same size OR one array and one scalar \n')
    return
end
%---------------------------------------------parser end -----------------------------------------------%
IE_suppression = 0.66019417475728148; % for inner electrodes (IE) voltage of 200V (or 100V?)
% CHANGE THAT

if isempty(AnchorBkg6G)
    AnchorBkg6G = GetAnchor6G('FPD_ROIlow',FPD_ROIlow);
end

%%switch Anchor6G
    %case 'OFF'
       % BKG = arrayfun(@(x,y) ((2314.07270639/((3.6./y)*x*1.e4)+72.87100681)/2.*2.86587134-131.37016528)/1000.*IE_suppression,MACE_Ba_T,WGTS_B_T);
    %case 'ON'
        BKG  = arrayfun(@(x,y) ((((2314.07270639/((3.6./y)*x*1.e4)+72.87100681)/2.*2.86587134-131.37016528)/1000.*IE_suppression )...
            ./ (((2314.07270639/((3.6./2.52)*6)+72.87100681)/2.*2.86587134-131.37016528)/1000.*IE_suppression) .* AnchorBkg6G),MACE_Ba_T,WGTS_B_T);
        %BKG  =  BKG ./148 .* numel(FPD_PixList);
%end
end


