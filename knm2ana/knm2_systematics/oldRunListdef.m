            % Run Summary version
            % Version = 'RunSummary-Prompt4b-fpd00';
            % Excluded Runs
            %             RunExcListIE       = unique([56175,56185,56283,56318,56408,56410,56631, 56687, 56705, 57021, 57037]);    % Inner electrode power supply breakdown
            %             RunExclListNaN     = [56331,56630];                                                                      % NaN value in K35
            %             RunExclListT2      = [56348,56632:56638,56665:56668,56675:56683];    % No T2 value or LARA problem
            %             RunExclListFPD     = [56280, 56304];             % detector slow control
            %             RunExclListqU      = [];[56160:56277,56332];        % voltage spike or other HV problems
            %             RunExclListE0      = 56343;                      % endpoint very low
            %             RunExclListOther   = [56415, 56604];             % other BREW comments
            %             RunExcList = [RunExcListIE,RunExclListNaN,RunExclListT2,RunExclListFPD,...
            %             RunExclListqU,RunExclListE0,RunExclListOther];
            
            
            %             % Read All KNM2 HD5 File
%             tmp = dir([getenv('SamakPath'), '/tritium-data/hdf5/',GetDataSet(FirstRun), '/*.h5']);
%             h5list = arrayfun(@(x) x.name,tmp,'UniformOutput',0);
%             h5list =  extractBefore(h5list,'.h5');
%             h5list = str2double(extractAfter(h5list,Version));
%             h5list(isnan(h5list)) = []; % %delete everything that has a different version
%             
%             % Truncate to exclude Runs from RunExcList:
%             h5list(h5list<FirstRun)=[];
%             h5list(h5list>LastRun)=[];
%             h5list(ismember(h5list,RunExcList))=[];