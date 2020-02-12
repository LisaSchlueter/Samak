function [PixList,RingPixList] = AziPatch2PixelCombi(myRingList,myPixList)
% function to synchronize list of pixel (myPixList) with list of azimuthal patches  
% input: - myPixList: list of  pixel, which you want to analyze
%        - myRingList:list of rings, which you want to analyze
% output: - PixList: list of selected pixels, which are within the selected rings
%         - RingPixList: cell of rings, which shows which pixel belongs to which rings

PixList_all = [];         
RingPixList = cell(5,1);  

% Center
if any(ismember(myRingList,1))
    POLE = [1:28];
    PixList_all = [PixList_all,POLE];
    RingPixList{1} = POLE;
end


if any(ismember(myRingList,2))
    NE = [28 29 30 40 41 42 52 53 54 64 65 66 76 77 78 88 89 90 100 101 102 112 113 114]+1;
    PixList_all = [PixList_all,NE];
    RingPixList{2} = NE;

end


if any(ismember(myRingList,3))
    SE = [37 38 39 49 50 51 61 62 63 73 74 75 85 86 87 97  98 99 109 110 111 121 122 123]+1;
    PixList_all = [PixList_all,SE];
    RingPixList{3} = SE;
end


if any(ismember(myRingList,4))
    SW = [34 35 36 46 47 48 58 59 60 70 71 72 82 83 84 94 95 96 106 107 108 118 119 120]+1;
    PixList_all = [PixList_all,SW];
    RingPixList{4} = SW;
end

if any(ismember(myRingList,5))
    NW = [31 32 33 43 44 45 55 56 57 67 68 69 79 80 81 91 92 93 103 104 106 115 116 117]+1;
    PixList_all = [PixList_all,NW];
    RingPixList{5} = NW;
end

% if any(ismember(myRingList,6))
%     OUT = [119:148];
%     PixList_all = [PixList_all,OUT];
%     RingPixList{6} = OUT ;
% end

 RingPixList = cellfun(@(x) intersect(myPixList,x),RingPixList,'UniformOutput',0); %cell, which shows which pixels belong to which ring
 PixList = intersect(myPixList,PixList_all);
 
end