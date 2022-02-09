

% actual FPD rings
[PixList,RingPixList] =Ring2Pixel(1:14,1:148);

p = PixList;
for i=1:numel(RingPixList)
    p(RingPixList{i}) = i;
end
FPDViewer(p)
%%

% 4 pseudo-rings
[PixList,RingPixList] =Ring2PixelCombi(1:14,1:148);
p = PixList;
for i=1:numel(RingPixList)
    p(RingPixList{i}) = i;
end
FPDViewer(p)