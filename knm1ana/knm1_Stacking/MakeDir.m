function MakeDir(dir)
if ~exist(dir,'dir')
    system(['mkdir -p ',dir]);
end
end