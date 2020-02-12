% convert into delta chi2 lookup table
dirName =  ([getenv('SamakPath'),'tritium-data/fit/TwinMCKnm1/']);
mydir = arrayfun(@(x) x.name,dir(dirName),'UniformOutput',0);
Index = cell2mat(cellfun(@(x) contains(x,'3000samples'),mydir,'UniformOutput',0));
%%
myfiles = mydir(Index);
d = cellfun(@(x) importdata(x),myfiles,'UniformOutput',0);

savedir = [getenv('SamakPath'),'tritium-data/FC/DeltaChi2LookupTable/']; 
system(['mkdir -p ',savedir]);
for i=1:numel(d)
   Chi2Best =  d{i}.Chi2Best;
   Chi2True = d{i}.Chi2True;
   DeltaChi2 = d{i}.DeltaChi2';
   mNuSqTrue = d{i}.mNuSq;
   mNuSqFit  = d{i}.par(1,:);
   
   filename = [savedir,strrep(myfiles{i},'TwinMC','DeltaChi2LookupTables_TwinMC')];
  
   save(filename,'Chi2Best','Chi2True','DeltaChi2','mNuSqTrue','mNuSqFit');

end