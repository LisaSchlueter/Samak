function  [plotHandle, cbHandle]  = PlotRingWiseFitPar(obj,varargin)
p=inputParser;
p.addParameter('PlotPar','Norm',@(x)ismember(x,{'qU','mTSq','Norm','Bkg'}));
p.addParameter('PlotParRef','',@(x)isfloat(x)); % plot parameter with respect to reference value (optional)
p.addParameter('SaveAs','',@(x)ischar(x) || isempty(x));
p.parse(varargin{:});
PlotPar = p.Results.PlotPar;
PlotParRef = p.Results.PlotParRef;
SaveAs = p.Results.SaveAs;

if ~strcmp(obj.AnaFlag,'Ring')
    fprintf('This routine makes only sense for multi-ring analysis \n')
    return
end

nRings = size(obj.RunData.qU,2);
switch PlotPar
    case 'qU'
        x = obj.FitResult.par(2*nRings+9:3*nRings+8);
        clabel = 'Retarding potential offset (eV)';
    case 'mTSq'
        x = obj.FitResult.par(3*nRings+10:4*nRings+9);
         clabel = sprintf('Energy smearing \\sigma^2 (eV^2)');
    case 'Norm'
        x = obj.FitResult.par(3+nRings:3+2*nRings-1);
         clabel = 'Normalization factor';
    case 'Bkg'
        x = obj.FitResult.par(3:2+nRings)*1e3;
         clabel = 'Background (mcps)';
end

if ~isempty(PlotParRef)
    x = x-PlotParRef;
end

xpixel = zeros(148,1).*NaN;
for i=1:nRings
    xpixel(obj.RingPixList{i}) = repmat(x(i),[numel(obj.RingPixList{i}),1]);
end

xpixel = xpixel-xpixel(1);
[plotHandle, cbHandle] = FPDViewer(xpixel);
cbHandle.Label.String = clabel;
cbHandle.Label.FontSize = get(gca,'FontSize')+4;


if ~isempty(SaveAs)
    export_fig(gca,'SaveAs');
end
end