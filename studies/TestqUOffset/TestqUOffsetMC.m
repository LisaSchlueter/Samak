mypath = [getenv('SamakPath'),'/studies/KNM1_RingAnalysis'];
range = 30;
if range==30
    excludedataStart = 17;
elseif range==90
    excludedataStart = 1;
end
RingList = 7:8;
nRings = numel(RingList);
FSDStart = 2.*nRings+3;
FSDStop  = 2*nRings+3+5;
qUStart  = 2*nRings+9;
qUStop   = 3*nRings+8;
fixPar   = ['',char(strjoin(string(FSDStart:FSDStop))),' ',char(strjoin(string(qUStart:qUStop)))];

%TwinBias_qU = 1:1:nRings;
TwinBias_Time = 60*60*24*1000;

CommonArg = {'RunNr',51443,'AnaFlag','Ring','RingList',RingList,'PullFlag',5,'fixPar',fixPar,...
             'DataType','Twin','exclDataStart',excludedataStart,'TwinBias_Time',TwinBias_Time,...
             'TwinBias_Bkg',[1 2]};%'TwinBias_qU',TwinBias_qU);

M = RunAnalysis(CommonArg{:});
%%
tic; M.Fit; toc;