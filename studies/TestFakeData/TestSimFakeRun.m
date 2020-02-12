clear all
nruns=360;

Runs{1}=SimFakeRun('RunNr',1,'RunTimeStart',datetime([2019 03 1 12 00 00]),'SimFakeRunPath','../../tritium-fakedata/mat/');
fprintf(2,'Fake Run %d Created \n ',1);

for run=2:1:nruns
    Runs{run-1}.RunTimeStart.Second = Runs{run-1}.RunTimeStart.Second + Runs{run-1}.TimeSec;
    Runs{run}=SimFakeRun('RunNr',run,'RunTimeStart',Runs{run-1}.RunTimeStart,'SimFakeRunPath','../../tritium-fakedata/mat/');
    fprintf(2,'Fake Run %d Created \n ',run);
end
variables    = fieldnames(Runs{1}); 
cell2mat(variables(2));

v=[];
for run=1:1:nruns
    v = [v Runs{run}.WGTS_MolFrac_DT];
end
%%
figure(1)
hfit = fitdist(v','Normal');disp(hfit);
histfit(v);PrettyFigureFormat;