
clear;
addpath(genpath('../../../Samak2.0'));

runn = 40531:41050;
kk = 1; jj = 1; ll = 1; nn = 1;
TotalTime = 0;
for ii = 1:length(runn)
    runname = ['2f-fpd00',num2str(runn(ii)),'.h5'];
    if exist(runname,'file')
        existingRun(kk) = runn(ii);
        kk = kk + 1;
        runData = load(num2str(runn(ii)));
        TimeX = runData.TimeSec;
        TotalTime = TotalTime + runData.TimeSec;
        if runData.TimeSec < 1*3600
            OneHoursRuns(nn,1) = runn(ii);
            nn = nn + 1;
        elseif runData.TimeSec < 1.8*3600
            OnePointFiveHooursScans(jj,1) = runn(ii);
            jj = jj + 1;  
        else
            ThreeHoursScans(ll,1) = runn(ii);
            ll = ll + 1;
        end
        
    end
end