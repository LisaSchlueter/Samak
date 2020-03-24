function myScanDir = GetKNM2ScanDirection(varargin)
% find scan direction of KNM2 tritium scans
% table from BREW (March 2020)
% 1 == up scan, 0 == down scan
p = inputParser;
p.addParameter('RunNr',56709,@(x)isfloat(x));
p.parse(varargin{:});
myRunNr = p.Results.RunNr;

filename = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/KNM2ScanDirection.txt'];
d = importdata(filename);

RunList   = d(:,1);
ScanDir = d(:,2);

myRunIndex = logical(sum(RunList==myRunNr,2));

if sum(myRunIndex)==0
    fprintf(2,'Invalid run number \n');
    return
end

myScanDir  = ScanDir(myRunIndex);
end