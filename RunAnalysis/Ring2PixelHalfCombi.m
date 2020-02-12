function [PixList,RingPixList] = Ring2PixelHalfCombi(myRingList,myPixList)
% function to synchronize list of pixel (myPixList) with list of rings (myRingList)
% input: - myPixList: list of  pixel, which you want to analyze
%        - myRingList:list of rings, which you want to analyze
% output: - PixList: list of selected pixels, which are within the selected rings
%          - RingPixList: cell of rings, which shows which pixel belongs to which rings

PixList_all = [];         
RingPixList = cell(2,1);  

if any(ismember(myRingList,1))
    PixList_all = [PixList_all,1:64];
    RingPixList{1}=1:64;
end

if any(ismember(myRingList,2))
    PixList_all = [PixList_all,65:148];
    RingPixList{2} =65:148;

end

 RingPixList = cellfun(@(x) intersect(myPixList,x),RingPixList,'UniformOutput',0); %cell, which shows which pixels belong to which ring
 PixList = intersect(myPixList,PixList_all);
 
end