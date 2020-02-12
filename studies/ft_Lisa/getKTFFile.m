function getKTFFile(runnr)
% RunSummaries
current_path = pwd;
RunSummary = sprintf('RunSummary2b-fpd00%u.ktf',runnr);
storeRS_here = sprintf('%s/../../tritium-data/ktf/RunSummaries/',current_path);

%PeriodSummaries
%PeriodSummary = 

cd /home/lisa/Dokumente/Studium/Masterstudium/Masterarbeit/RS2HDF5/install/bin;
command = sprintf('idle GetFile %s@FirstTritium.katrin > %s/%s',RunSummary,storeRS_here, RunSummary);
system(command);
end