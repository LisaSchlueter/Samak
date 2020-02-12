% -----------------------------------------------------------------------------------------------
% Function, which returns a RunList (e.g. for MultiRunAnalysis)
% input: Stackfilenae:
%       -'StackCD100all'    : all first tritium scans with 100% column density
%       -'StackCD100up '    : all first tritium UP scans with 100% column density
%       -'StackCD100down'   : all first tritium DOWN scans with 100% column density
%       -'StackCD100random' : all first tritium RANDOM scans with 100% column density
% -----------------------------------------------------------------------------------------------
% Lisa Schlüter
% TUM/MPP
% 09/2018
% -----------------------------------------------------------------------------------------------
function [RunList] = GetFTRunList(StackFileName)     
% First tritium runs lists 
            %(taken from BREW, http://katana.npl.washington.edu/~sanshiro/brew/FirstTritiumTest/)
            RunList_100up = [...% 100% column density - up scans (== from -1600eV to -1800eV)
                40539, 40541, 40543, 40604, 40611, 40613, 40667, 40669, 40671, 40673, ...
                40675, 40677, 40679, 40681, 40683, 40685, 40687, 40689, 40691, 40693, ...
                40977, 40980, 40983, 40986, 40989, 40992, 40995, 41002, 41005, 41008, ...
                41011, 41014, 41019, 41022,  41026,  41029, 41032];
            RunList_100down = [...% 100% column density - down scans (== from -1800eV to -1600eV)
                40531, 40540, 40542, 40603, 40610, 40612, 40668, 40670, 40672, 40674, ...
                40676, 40678, 40680, 40682, 40684, 40686, 40688, 40690, 40692, 40976, ...
                40979, 40982, 40985, 40988, 40991, 40994, 40997, 41007, 41010, 41013, ...
                41016, 41017, 41020, 41023, 41025,  41028, 41031];
            RunList_100random = [...% 100% column density - random scans
                41003,41006, 41015, 41021, 41024,  41027,  41030, 41033];
            switch StackFileName
                case 'StackCD100all'%all runs (up+down scans) with 100% column density
                    RunList = sort([RunList_100up,RunList_100down]);
                case 'StackCD100up'
                    RunList = RunList_100up;
                case 'StackCD100down'
                    RunList = RunList_100down;
                case 'StackCD100random'
                    RunList = RunList_100random;
                case 'StackCD100_3hours'  
                    RunList_CD100all = sort([RunList_100up,RunList_100down]);
                    start3h = find(RunList_CD100all==40667);
                    stop3h  = find(RunList_CD100all==40693);
                    RunList = RunList_CD100all(start3h:stop3h);
                otherwise
                    fprintf('RunList Name unknown! \n');
            end
end