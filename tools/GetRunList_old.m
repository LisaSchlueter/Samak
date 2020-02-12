% function [ n, datafiles ] = GetRunList( directory, filepattern, filterout, fpattern)
%
% input:
% - directory
% - filepattern
% - filterout | 1 or 0
% - fpattern | not include files with this pattern
% fpattern is charracter array created by fpattern =
% char('string1', 'string2', ...,'stringn');
% or in case of one string fpattern = 'string';
%
% output:
% - n | number of files
% - datafiles | structure with datafiles
% datafiles(i).name -- sigle file name
% datafiles(i).fullname -- including full path
% datafiles(i).runlist -- runlist

function [ n, datafiles, runlist] = GetRunList( directory, filepattern, filterout, fpattern)
%initialisation of the structure
datafiles(1).name     = '';
datafiles(1).fullname = '';
datafiles(1).runlist  = '';
runlist=[];
if nargin <3
    filterout = 0;
end
tmpdir = dir([directory filepattern '*']);
nn = length(tmpdir);
%filtering not necessary files
cnt = 0;
disp('Folder contain files for further processing:')
if filterout==1
    [nstr s] = size(fpattern);
    filter_flag = 0;
    for i=1:nn
        for st=1:nstr
            S = regexp(tmpdir(i).name,fpattern(st,:));
            if(isempty(S))
                filter_flag=0;
            else
                filter_flag=1;
            end
        end
        if(filter_flag==0)
            % filter pattern not maching ... so this is correct file
            cnt = cnt+1;
            datafiles(cnt).name     = tmpdir(i).name;
            datafiles(cnt).fullname = [directory tmpdir(i).name];
            [~,run,~] = fileparts(datafiles(cnt).name); datafiles(cnt).runlist  = run;
            %run = extractAfter(run,"2f-fpd00");
            run = extractBefore(run,"ex2");
            tmpdir(i).name;
            disp(run);
            runlist = [runlist str2num(run)];
        else
            %filter pattern matching
            %display(['Filtering out: ' tmpdir(i).name ' file'])
        end
    end
else
    for i=1:nn
        datafiles(i).name     = tmpdir(i).name;
        datafiles(i).fullname = [directory tmpdir(i).name];
        [~,run,~] = fileparts(datafiles(cnt).name); datafiles(cnt).runlist  = run;
        runlist = [runlist run];
    end
end
n = length(datafiles);
for i=1:n
    display([ num2str(i) ' file: ' datafiles(i).name ' run: ' datafiles(i).runlist])
end
display(['Number of listed files: ' num2str(n)]);
return