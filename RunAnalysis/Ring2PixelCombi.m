function [PixList,RingPixList] = Ring2PixelCombi(myRingList,myPixList)
% function to synchronize list of pixel (myPixList) with list of rings (myRingList)
% input: - myPixList: list of  pixel, which you want to analyze
%        - myRingList:list of rings, which you want to analyze
% output: - PixList: list of selected pixels, which are within the selected rings
%          - RingPixList: cell of rings, which shows which pixel belongs to which rings

PixList_all = [];         
RingPixList = cell(4,1);  

if any(ismember(myRingList,1))
    PixList_all = [PixList_all,1:28];
    RingPixList{1}=1:28;
end


if any(ismember(myRingList,2))
    PixList_all = [PixList_all,29:64];
    RingPixList{2} =29:64;
 %   RingPixList{2} =29:40;

end


if any(ismember(myRingList,3))
    PixList_all = [PixList_all,65:100];
    RingPixList{3} =65:100;
end


if any(ismember(myRingList,4))
    PixList_all = [PixList_all,101:148];
    RingPixList{4} =101:148;
   % RingPixList{4} =101:112;

end


 RingPixList = cellfun(@(x) intersect(myPixList,x),RingPixList,'UniformOutput',0); %cell, which shows which pixels belong to which ring
 PixList = intersect(myPixList,PixList_all);
 
end