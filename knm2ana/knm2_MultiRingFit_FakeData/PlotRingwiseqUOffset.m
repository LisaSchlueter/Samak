function  [plotHandle, cbHandle]  = PlotRingwiseqUOffset(obj,varargin)
subRun = 1;
qUpixel = zeros(148,1).*NaN;

switch obj.AnaFlag
    case 'Ring'
 for i=1:size(obj.RunData.qU,2)
     qUpixel(obj.RingPixList{i}) = repmat(obj.RunData.qU(subRun,i),[numel(obj.RingPixList{i}),1]);
 end
    case 'StackPixel'
     qUpixel(obj.PixList) = repmat(obj.RunData.qU(subRun),[numel(obj.PixList),1]); 
end
 
 qUpixel = qUpixel-qUpixel(1);
 [plotHandle, cbHandle] = FPDViewer(qUpixel);
 cbHandle.Label.String = 'Retarding potential offset (eV)';
 cbHandle.Label.FontSize = get(gca,'FontSize')+4;
end