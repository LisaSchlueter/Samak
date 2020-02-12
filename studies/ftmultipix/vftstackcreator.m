%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endpoint sensitivity study
%
% Configuration for First Tritium May

% P. Morales 2018
% Last update 23/04/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%Initialization
addpath(genpath('../../../Samak2.0'));

DTconcen = 'all'; mpix_suffix = 'mpix';

if strcmp(mpix_suffix,'mpix')
    r = (1:8)';
    x = repmat((1:148)',1,8);
else
    r = ones(8,1);
    x = (1:8);
end

switch DTconcen
    case '9'
        stackfilename1 = ['40257',mpix_suffix,'.mat'];
        stackfilename2 = ['40258',mpix_suffix,'.mat'];
        stackfilename3 = ['40259',mpix_suffix,'.mat'];
        stackfilename4 = ['40260',mpix_suffix,'.mat'];
    case '12'
        stackfilename1 = ['40263',mpix_suffix,'.mat'];
        stackfilename2 = ['40264',mpix_suffix,'.mat'];
        stackfilename3 = ['40265',mpix_suffix,'.mat'];
        stackfilename4 = ['40266',mpix_suffix,'.mat'];
    case 'all'
        stackfilename1 = ['40257',mpix_suffix,'.mat'];
        stackfilename2 = ['40258',mpix_suffix,'.mat'];
        stackfilename3 = ['40259',mpix_suffix,'.mat'];
        stackfilename4 = ['40260',mpix_suffix,'.mat'];
        stackfilename5 = ['40263',mpix_suffix,'.mat'];
        stackfilename6 = ['40264',mpix_suffix,'.mat'];
        stackfilename7 = ['40265',mpix_suffix,'.mat'];
        stackfilename8 = ['40266',mpix_suffix,'.mat'];
end

R1 = load(stackfilename1);
R2 = load(stackfilename2);
R3 = load(stackfilename3);
R4 = load(stackfilename4);

if strcmp(DTconcen,'all')
    R5 = load(stackfilename5);
    R6 = load(stackfilename6);
    R7 = load(stackfilename7);
    R8 = load(stackfilename8);
end



if size(R1.qU,1) == 28
    R1.qU = R1.qU(3:28,:);
    R1.qUfrac = R1.qUfrac(3:28,:);
    R1.TBDIS = R1.TBDIS(3:28,:);
end

if size(R2.qU,1) == 28
    R2.qU = R2.qU(3:28,:);
    R2.qUfrac = R2.qUfrac(3:28,:);
    R2.TBDIS = R2.TBDIS(3:28,:);
end

msize = size(R1.qU,2);

limit = 110;

% Stack qU
qU57 = reshape(R1.qU(R1.qU <= (18575 - limit)),[12 msize]);

qUall_M(:,x(:,1),r(1)) = reshape(R1.qU(R1.qU > (18575 - limit)),[14 msize]);
qUall_M(:,x(:,2),r(2)) = reshape(R2.qU(R2.qU > (18575 - limit)),[14 msize]);
qUall_M(:,x(:,3),r(3)) = reshape(R3.qU(R3.qU > (18575 - limit)),[14 msize]);
qUall_M(:,x(:,4),r(4)) = reshape(R4.qU(R4.qU > (18575 - limit)),[14 msize]);

if strcmp(DTconcen,'all')
    qUall_M(:,x(:,5),r(5)) = reshape(R5.qU(R5.qU > (18575 - limit)),[14 msize]);
    qUall_M(:,x(:,6),r(6)) = reshape(R6.qU(R6.qU > (18575 - limit)),[14 msize]);
    qUall_M(:,x(:,7),r(7)) = reshape(R7.qU(R7.qU > (18575 - limit)),[14 msize]);
    qUall_M(:,x(:,8),r(8)) = reshape(R8.qU(R8.qU > (18575 - limit)),[14 msize]);
end
if strcmp(mpix_suffix,'mpix')
    qUall = mean(qUall_M,3);
else
    qUall = mean(qUall_M,2);
end

qU = [qU57;qUall];

% stack TBDIS
TBD57 = reshape(R1.TBDIS(R1.qU <= (18575 - 110)),[12 msize]);

TBDall_M(:,x(:,1),r(1)) = reshape(R1.TBDIS(R1.qU > (18575 - limit)),[14 msize]);
TBDall_M(:,x(:,2),r(2)) = reshape(R2.TBDIS(R2.qU > (18575 - limit)),[14 msize]);
TBDall_M(:,x(:,3),r(3)) = reshape(R3.TBDIS(R3.qU > (18575 - limit)),[14 msize]);
TBDall_M(:,x(:,4),r(4)) = reshape(R4.TBDIS(R4.qU > (18575 - limit)),[14 msize]);

if strcmp(DTconcen,'all')
    TBDall_M(:,x(:,5),r(5)) = reshape(R5.TBDIS(R5.qU > (18575 - limit)),[14 msize]);
    TBDall_M(:,x(:,6),r(6)) = reshape(R6.TBDIS(R6.qU > (18575 - limit)),[14 msize]);
    TBDall_M(:,x(:,7),r(7)) = reshape(R7.TBDIS(R7.qU > (18575 - limit)),[14 msize]);
    TBDall_M(:,x(:,8),r(8)) = reshape(R8.TBDIS(R8.qU > (18575 - limit)),[14 msize]);
end

if strcmp(mpix_suffix,'mpix')
    TBDall = sum(TBDall_M,3);
else
    TBDall = sum(TBDall_M,2);
end

TBDIS = [TBD57;TBDall];

% stack runtimes

runtime1 = R1.TimeSec*R1.qUfrac;
runtime2 = R2.TimeSec*R2.qUfrac;
runtime3 = R3.TimeSec*R3.qUfrac;
runtime4 = R4.TimeSec*R4.qUfrac;

if strcmp(DTconcen,'all')
    runtime5 = R5.TimeSec*R6.qUfrac;
    runtime6 = R6.TimeSec*R7.qUfrac;
    runtime7 = R7.TimeSec*R8.qUfrac;
    runtime8 = R8.TimeSec*R8.qUfrac;
end

runtime57 = runtime1(R1.qU(:,1) <= (18575 - limit));

runtime_M(:,1) = runtime1(R1.qU(:,1) > (18575 - limit));
runtime_M(:,2) = runtime2(R2.qU(:,1) > (18575 - limit));
runtime_M(:,3) = runtime3(R3.qU(:,1) > (18575 - limit));
runtime_M(:,4) = runtime4(R4.qU(:,1) > (18575 - limit));

if strcmp(DTconcen,'all')
    runtime_M(:,5) = runtime5(R5.qU(:,1) > (18575 - limit));
    runtime_M(:,6) = runtime6(R6.qU(:,1) > (18575 - limit));
    runtime_M(:,7) = runtime7(R7.qU(:,1) > (18575 - limit));
    runtime_M(:,8) = runtime8(R8.qU(:,1) > (18575 - limit));
end

runtimeall = sum(runtime_M,2);
stackruntime = [runtime57;runtimeall];
TimeSec = sum(stackruntime);
qUfrac = stackruntime/sum(stackruntime);

% stack DT frac

DTfrac(1) = R1.WGTS_MolFrac_DT;
DTfrac(2) = R2.WGTS_MolFrac_DT;
DTfrac(3) = R3.WGTS_MolFrac_DT;
DTfrac(4) = R4.WGTS_MolFrac_DT;

if strcmp(DTconcen,'all')
    DTfrac(5) = R5.WGTS_MolFrac_DT;
    DTfrac(6) = R6.WGTS_MolFrac_DT;
    DTfrac(7) = R7.WGTS_MolFrac_DT;
    DTfrac(8) = R8.WGTS_MolFrac_DT;
end

WGTS_MolFrac_DT = mean(DTfrac);
WGTS_MolFrac_TT = 0;
WGTS_MolFrac_HT = 0;

% column density stack

CD(1) = R1.WGTS_CD_MolPerCm2;
CD(2) = R2.WGTS_CD_MolPerCm2;
CD(3) = R3.WGTS_CD_MolPerCm2;
CD(4) = R4.WGTS_CD_MolPerCm2;

if strcmp(DTconcen,'all')
    CD(5) = R5.WGTS_CD_MolPerCm2;
    CD(6) = R6.WGTS_CD_MolPerCm2;
    CD(7) = R7.WGTS_CD_MolPerCm2;
    CD(8) = R8.WGTS_CD_MolPerCm2;
end

WGTS_CD_MolPerCm2 = mean(CD);

% stack background

if strcmp(mpix_suffix,'mpix')
    TBDIS_bck = sum(TBDIS((end-2:end),:));
    TBDIS(end-2,:) = TBDIS_bck;
    TBDIS((end-1:end),:) = [];
    qUfrac(end-2) = sum(qUfrac(end-2:end));
    qUfrac(end-1:end) = [];
    qU(end-1:end,:) = [];
end




% save(['../../inputs/ftdata/stackDT',DTconcen,'.mat'],...
%     'stackTBDIS','stackqU','stackTimeSec','stackqUfrac','WGTS_CD_MolPerCm2',...
%     'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
%     '-v7.3','-nocompression')

save(['../../tritium-data/stackDT',DTconcen,mpix_suffix,'.mat'],...
    'TBDIS','qU','TimeSec','qUfrac','WGTS_MolFrac_DT',...
    'WGTS_MolFrac_TT','WGTS_MolFrac_HT',...
    'TimeSec','WGTS_CD_MolPerCm2','-v7.3','-nocompression')

RunTime = TimeSec;
save(['../../simulation/katrinsetup/TD_DataBank/stackRunDT',...
    DTconcen,mpix_suffix,'.mat'],...
    'qU','qUfrac','RunTime','-v7.3','-nocompression')


