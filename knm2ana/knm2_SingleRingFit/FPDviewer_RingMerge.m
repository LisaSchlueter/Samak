
DataSet = 'Knm2';
RingMerge = 'None';

PixList = GetPixList(DataSet);
RingList = 1:12; %exclude last ring

switch RingMerge
    case 'Default'
        [PixList,RingPixList] = Ring2PixelDefCombi(RingList,PixList);
        RingList = 1:10;
    case 'None'
        [PixList,RingPixList] = Ring2Pixel(RingList,PixList);
        RingList = 1:12;
    case 'Full'
        [PixList,RingPixList] = Ring2PixelCombi(RingList,PixList);
        RingList = 1:4;%1:numel();
    case 'Half'
        [PixList,RingPixList] = Ring2PixelHalfCombi(RingList,PixList);
        RingList = 1:2;
    case 'Azi'
        [PixList,RingPixList] = AziPatch2PixelCombi(RingList,PixList);
        RingList = 1:5;
    case {'AziHalfNS','AziHalfEW'}
        [PixList,RingPixList] = AziHalfPatch2PixelCombi(RingList,PixList,RingMerge);
        RingList = 1:2;
    case {'Slice','Slice2','Slice3','Slice4','Slice3_1'}
        [PixList,RingPixList] = Slice2PixelCombi(PixList,RingMerge);
        RingList = 1:numel(RingPixList); % now psuedo-rings (before, no meaning in slice mode)
end

FPD_pix = NaN.*zeros(148,1);
for i=1:numel(RingList)
    FPD_pix(RingPixList{i}) = i;
end


[~, cbHandle,~] = FPDViewer(FPD_pix);
cbHandle.delete;
colormap(jet);
xlabel(RingMerge)