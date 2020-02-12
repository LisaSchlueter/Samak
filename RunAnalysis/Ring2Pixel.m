function [PixList,RingPixList] = Ring2Pixel(myRingList,myPixList)
% function to synchronize list of pixel (myPixList) with list of rings (myRingList)
% input: - myPixList: list of  pixel, which you want to analyze
%        - myRingList:list of rings, which you want to analyze
% output: - PixList: list of selected pixels, which are within the selected rings
%          - RingPixList: cell of rings, which shows which pixel belongs to which rings

PixList_all = [];         
RingPixList = cell(12,1);  

if any(ismember(myRingList,1))
    PixList_all = [PixList_all,1:4];
    RingPixList{1}=1:4;
end

if any(ismember(myRingList,2))
    PixList_all = [PixList_all,5:16];
    RingPixList{2} =5:16;
end

if any(ismember(myRingList,3))
    PixList_all = [PixList_all,17:28];
    RingPixList{3} =17:28;
end

if any(ismember(myRingList,4))
    PixList_all = [PixList_all,29:40];
    RingPixList{4} =29:40;
end

if any(ismember(myRingList,5))
    PixList_all = [PixList_all,41:52];
    RingPixList{5} =41:52;
end

if any(ismember(myRingList,6))
    PixList_all = [PixList_all,53:64];
    RingPixList{6} =53:64;
end

if any(ismember(myRingList,7))
    PixList_all = [PixList_all,65:76];
    RingPixList{7} =65:76;
end

if any(ismember(myRingList,8))
    PixList_all = [PixList_all,77:88];
    RingPixList{8} =77:88;
end

if any(ismember(myRingList,9))
    PixList_all = [PixList_all,89:100];
    RingPixList{9} =89:100;
end

if any(ismember(myRingList,10))
    PixList_all = [PixList_all,101:112];
    RingPixList{10} =101:112;
end

if any(ismember(myRingList,11))
    PixList_all = [PixList_all,113:124];
    RingPixList{11} =113:124;
end

if any(ismember(myRingList,12))
    PixList_all = [PixList_all,125:136];
    RingPixList{12} =125:136;
end

if any(ismember(myRingList,13))
    PixList_all = [PixList_all,137:148];
    RingPixList{13} =137:148;
end

 RingPixList = cellfun(@(x) intersect(myPixList,x),RingPixList,'UniformOutput',0); %cell, which shows which pixels belong to which ring
 PixList = intersect(myPixList,PixList_all);
end