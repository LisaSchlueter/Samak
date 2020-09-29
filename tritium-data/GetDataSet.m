function DataSet = GetDataSet(RunNr)
if isempty(RunNr)
    DataSet = '';
    % do nothing
elseif all(~isletter(RunNr)) % Single Run
    
    if ischar(RunNr)
        RunNr = str2double(RunNr);
    end
    
    if all(RunNr <= 1000)
        DataSet = 'Knm1';
    elseif all(RunNr <= 10000)
        DataSet = 'Knm2';
    elseif all(RunNr <= 41033)
        DataSet = 'FirstTritium.katrin';
    elseif all(RunNr <= 51936)
        DataSet = 'Knm1';
    elseif all(RunNr <= 57136)
        DataSet = 'Knm2';
    elseif all(RunNr <= 62957)
        DataSet = 'Knm3a';
    else
        DataSet = 'Knm3b';
    end
    
    %MultiRuns
elseif contains(RunNr,'FirstTritium') || ismember(RunNr,{'StackCD100all','StackCD100_3hours','FTpaper'})
    DataSet = 'FirstTritium.katrin';
elseif contains(RunNr,'KNM1')
    DataSet = 'Knm1';
elseif contains(RunNr,'KNM2')
    DataSet = 'Knm2';
elseif contains(RunNr,'Stack') %for generic labels
    RunNr = strrep(RunNr,'Stack','');
    RunNr = str2double(strtok(RunNr,'_')); % Get first run
    if RunNr <= 41033
        DataSet = 'FirstTritium.katrin';
    elseif all(RunNr <= 51936)
        DataSet = 'Knm1';
    else
        DataSet = 'Knm2';
    end
end
end