function FPD_Coverage = GetFPDcoverage(varargin)
%-----------------parser start ---------------------------------
p=inputParser;
p.addParameter('FPD_ROIlow',26,@(x)isfloat(x)); % Region of Interest (ROI)- lower cut [keV]
p.addParameter('FPD_ROIup',32,@(x)isfloat(x));  % fixed and not used at the moment
p.parse(varargin{:});
FPD_ROIlow = p.Results.FPD_ROIlow;
FPD_ROIup  = p.Results.FPD_ROIup;
%-----------------parser end ---------------------------------

% FPD Coverage versus ROI - Sanshiro data, 31/10/2018
FPDfile = load('FPDROICoverage.txt');  %FPD covarage as a function of ROIlow
FPD_Coverage = interp1(FPDfile(:,1),FPDfile(:,2),FPD_ROIlow);
end