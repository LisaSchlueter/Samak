function Anchor6G = GetAnchor6G(varargin)
%-----------------parser start ---------------------------------
p=inputParser;
p.addParameter('FPD_ROIlow',26,@(x)isfloat(x)); % Region of Interest (ROI)- lower cut [keV]
p.addParameter('FPD_ROIup',32,@(x)isfloat(x));  % fixed and not used at the moment
p.parse(varargin{:});
FPD_ROIlow = p.Results.FPD_ROIlow;
FPD_ROIup  = p.Results.FPD_ROIup;
%-----------------parser end ---------------------------------

% Background coverage versus ROI -  Anna data, 31/10/2018
BKGfile = load('BackgroundROICut.txt');
Anchor6G = interp1(BKGfile(:,1),BKGfile(:,2)*1000,FPD_ROIlow)*1e-03; %Background as a function of ROIlow
end