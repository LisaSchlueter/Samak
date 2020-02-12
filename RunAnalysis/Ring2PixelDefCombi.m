function [PixList,RingPixList] = Ring2PixelDefCombi(myRingList,myPixList)
PixList_all = [];         
RingPixList = cell(10,1);  


if any(ismember(myRingList,1))
    PixList_all = [PixList_all,1:16];
    RingPixList{1} =1:16;
end

if any(ismember(myRingList,2))
    PixList_all = [PixList_all,17:28];
    RingPixList{2} =17:28;
end

if any(ismember(myRingList,3))
    PixList_all = [PixList_all,29:40];
    RingPixList{3} =29:40;
end

if any(ismember(myRingList,4))
    PixList_all = [PixList_all,41:52];
    RingPixList{4} =41:52;
end

if any(ismember(myRingList,5))
    PixList_all = [PixList_all,53:64];
    RingPixList{5} =53:64;
end

if any(ismember(myRingList,6))
    PixList_all = [PixList_all,65:76];
    RingPixList{6} =65:76;
end

if any(ismember(myRingList,7))
    PixList_all = [PixList_all,77:88];
    RingPixList{7} =77:88;
end

if any(ismember(myRingList,8))
    PixList_all = [PixList_all,89:100];
    RingPixList{8} =89:100;
end

if any(ismember(myRingList,9))
    PixList_all = [PixList_all,101:112];
    RingPixList{9} =101:112;
end

if any(ismember(myRingList,10))
    PixList_all = [PixList_all,113:148];
    RingPixList{10} =113:148;
end

 RingPixList = cellfun(@(x) intersect(myPixList,x),RingPixList,'UniformOutput',0); %cell, which shows which pixels belong to which ring
 PixList = intersect(myPixList,PixList_all);
 
end