function StackRuns(obj,varargin)

% Init
obj.RunsData = cell(obj.nRuns,1);

for rr = 1:obj.nRuns
    obj.RunsData{rr}   = load([num2str(obj.RunList(rr)),obj.mpix,obj.ringCutFlag,'.mat']);
    RunsqU(:,:,rr) = obj.RunsData{rr}.qU; %#ok<AGROW>
end

qUmedian = median(RunsqU,3);

StackCounts = zeros(size(qUmedian));       % Sum of Counts per Subrun
StackqU = zeros(size(qUmedian));           % average retarding potential
StackqUFrac  = zeros(size(qUmedian,1),1);  % average time per Subrun
StackTime = 0;                             % total time
StackCD = 0;                               % average column density
StackCDSubRun = zeros(size(qUmedian,1),1); % average column density per Subrun
StackTT = 0; StackDT = 0; StackHT = 0;     % average Isotopologue fractions
StackDTSubRun = zeros(size(qUmedian,1),1); % average fractions per Subrun
StackHTSubRun = zeros(size(qUmedian,1),1);
StackTTSubRun = zeros(size(qUmedian,1),1);
obj.StackFileName = 'Stack_FT';

% First Iteration: rough selection
nsr = 1;
for rr = 1:obj.nRuns
    if ~all(all(abs(obj.RunsData{rr}.qU - qUmedian) < obj.StackTolerance))
        obj.NotStackedRuns(nsr) = obj.RunList(rr);
        nsr = nsr + 1;
    end
end
qUmedian = median(RunsqU(:,:,~ismember(obj.RunList,obj.NotStackedRuns)),3);
% Second Iteration: fine selection
sr = 1;
for rr = 1:obj.nRuns
    if all(all(abs(obj.RunsData{rr}.qU - qUmedian) < obj.StackTolerance))
        obj.StackedRuns(sr) = obj.RunList(rr);
        % Sum Counts and Time
%         if all(~ismember(obj.DataEffCorr,{'OFF'}))
%             obj.ROIPileUpDataCorrectionStackRuns(rr);
%         end
        % Monte Carlo data per Run
        Model_MC = ref_RunAnalysis(obj.RunList(rr),obj.ringCutFlag,obj.mpix,'ISCS','Theory','recomputeRF','OFF');
        Model_MC.ComputeTBDDS; Model_MC.ComputeTBDIS;
        Model_MC.AddStatFluctTBDIS();
        %DatatoSave = [Model_MC.qU,Model_MC.TBDIS];
        %save(['MCRunData/',num2str(obj.RunList(rr)),'.mat'],DatatoSave);
        obj.RunsData{rr}.TBDIS = Model_MC.TBDIS;
        StackCounts   = StackCounts + obj.RunsData{rr}.TBDIS;
        StackTime     = StackTime + obj.RunsData{rr}.TimeSec;
        % Average everything else
        StackqU       = StackqU + obj.RunsData{rr}.qU*obj.RunsData{rr}.TimeSec;
        StackqUFrac   = StackqUFrac + obj.RunsData{rr}.qUfrac*obj.RunsData{rr}.TimeSec;
        StackCD       = StackCD + obj.RunsData{rr}.WGTS_CD_MolPerCm2*obj.RunsData{rr}.TimeSec;
        StackCDSubRun = StackCDSubRun + obj.RunsData{rr}.WGTS_CD_MolPerCm2_SubRun*obj.RunsData{rr}.TimeSec;
        StackTT       = StackTT + obj.RunsData{rr}.WGTS_MolFrac_TT*obj.RunsData{rr}.TimeSec;
        StackDT       = StackDT + obj.RunsData{rr}.WGTS_MolFrac_DT*obj.RunsData{rr}.TimeSec;
        StackHT       = StackHT + obj.RunsData{rr}.WGTS_MolFrac_HT*obj.RunsData{rr}.TimeSec;
        StackDTSubRun = StackDTSubRun + obj.RunsData{rr}.WGTS_MolFrac_DT_SubRun*obj.RunsData{rr}.TimeSec;
        StackHTSubRun = StackHTSubRun + obj.RunsData{rr}.WGTS_MolFrac_HT_SubRun*obj.RunsData{rr}.TimeSec;
        StackTTSubRun = StackTTSubRun + obj.RunsData{rr}.WGTS_MolFrac_TT_SubRun*obj.RunsData{rr}.TimeSec;
        % Labeling
        %RunString           = num2str(obj.RunList(rr));
        %obj.StackFileName   = [obj.StackFileName,'_',RunString(3:end)];
        
        sr = sr + 1;
        
    else
        obj.NotStackedRuns(nsr) = obj.RunList(rr);
        nsr = nsr + 1;
    end
end

% Prepare variable to save data
TBDIS                    = StackCounts; %#ok
TimeSec                  = StackTime;
qU                       = StackqU./TimeSec; %#ok
qUfrac                   = StackqUFrac./TimeSec; %#ok

WGTS_CD_MolPerCm2        = StackCD/TimeSec; %#ok
WGTS_MolFrac_TT          = StackTT/TimeSec; %#ok
WGTS_MolFrac_DT          = StackDT/TimeSec; %#ok
WGTS_MolFrac_HT          = StackHT/TimeSec; %#ok
WGTS_CD_MolPerCm2_SubRun = StackCDSubRun./TimeSec; %#ok
WGTS_MolFrac_DT_SubRun   = StackDTSubRun./TimeSec; %#ok
WGTS_MolFrac_HT_SubRun   = StackHTSubRun./TimeSec; %#ok
WGTS_MolFrac_TT_SubRun   = StackTTSubRun./TimeSec; %#ok

save(['fake-dataFT/',obj.StackFileName,obj.mpix,obj.ringCutFlag,'.mat'],...
    'TBDIS','qU','TimeSec','qUfrac',...
    'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
    'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
    'WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun','WGTS_MolFrac_TT_SubRun',...
    '-v7.3','-nocompression')

obj.RunData = load([obj.StackFileName,obj.mpix,obj.ringCutFlag,'.mat']);
fprintf(2,'Showing Stacked Data (data for Ring selected later) \n')
disp(obj.RunData);

% Save stacked pixels excluding two outer rings
qU = mean(qU(:,1:124),2); 
TBDIS = sum(TBDIS(:,1:124),2);

save(['fake-dataFT/',obj.StackFileName,'ex2.mat'],...
    'TBDIS','qU','TimeSec','qUfrac',...
    'WGTS_CD_MolPerCm2','WGTS_CD_MolPerCm2_SubRun',...
    'WGTS_MolFrac_TT','WGTS_MolFrac_DT','WGTS_MolFrac_HT',...
    'WGTS_MolFrac_TT_SubRun','WGTS_MolFrac_DT_SubRun','WGTS_MolFrac_HT_SubRun',...
    '-v7.3','-nocompression')


% switch obj.AnaFlag
%     case 'MultiPixel'
%         obj.RunData.TBDIS = obj.RunData.TBDIS(:,obj.PixList);
%         obj.RunData.qU = obj.RunData.qU(:,obj.PixList);
%     case 'SinglePixel'
%         obj.RunData.TBDIS = obj.RunData.TBDIS(:,obj.PixList(1));
%         obj.RunData.qU = obj.RunData.qU(:,obj.PixList(1));
%     case 'Ring'
%         % Data is selected in Simulate Model, where the Model
%         % has the information on the relation pixel-ring
%         
% end

end % StackRuns
