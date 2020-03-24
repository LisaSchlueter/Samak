function [PixList,RingPixList] = AziHalfPatch2PixelCombi(myRingList,myPixList,RingMerge)
% function to synchronize list of pixel (myPixList) with list of azimuthal patches  
% input: - myPixList: list of  pixel, which you want to analyze
%        - myRingList:list of rings, which you want to analyze
% output: - PixList: list of selected pixels, which are within the selected rings
%         - RingPixList: cell of rings, which shows which pixel belongs to which rings

PixList_all = [];
RingPixList = cell(2,1);

% define pixel lists
PoleNE = [0,5:7,16:18]+1;                                                               % pole north east
NE = [28 29 30 40 41 42 52 53 54 64 65 66 76 77 78 88 89 90 100 101 102 112 113 114]+1; % north east
PoleSE = [3,4,14,15,25:27]+1;                                                           % pole south east
SE = [37 38 39 49 50 51 61 62 63 73 74 75 85 86 87 97  98 99 109 110 111 121 122 123]+1;% south east
PoleNW = [1,8:10,  19:21]+1;                                                            % pole north west
PoleSW = [2,11:13, 22:24]+1;                                                            % pole south west
NW = [31 32 33 43 44 45 55 56 57 67 68 69 79 80 81 91 92 93 103 104 106 115 116 117]+1; % north west
SW = [34 35 36 46 47 48 58 59 60 70 71 72 82 83 84 94 95 96 106 107 108 118 119 120]+1; % south west

% ring 1: first half
if any(ismember(myRingList,1))
    switch RingMerge
        case 'AziHalfNS' % ring 1: north
            PixList_all = [PixList_all,NE,NW,PoleNE,PoleNW];
            RingPixList{1} = [NE,NW,PoleNE,PoleNW];
        case 'AziHalfEW' % ring 1: east
            PixList_all = [PixList_all,NE,SE,PoleNE,PoleSE];
            RingPixList{1} = [NE,SE,PoleNE,PoleSE];
    end
end

% ring 2: second half
if any(ismember(myRingList,2))
    switch RingMerge
        case 'AziHalfNS' % ring 2: south
            PixList_all = [PixList_all,SE,SW,PoleSE,PoleSW];
            RingPixList{2} = [SE,SW,PoleSE,PoleSW];
        case 'AziHalfEW' % ring 1: west
            PixList_all = [PixList_all,NW,SW,PoleNW,PoleSW];
            RingPixList{2} = [NW,SW,PoleNW,PoleSW];
    end
end

 RingPixList = cellfun(@(x) intersect(myPixList,x),RingPixList,'UniformOutput',0); %cell, which shows which pixels belong to which ring
 PixList = intersect(myPixList,PixList_all);
end