% little script to add new RS parameters to old MC twins
% 1. Bmax
% 2. RW_BiasVoltage

RunAnaArg = {'RunList','KNM2_RW3',...  % define run number -> see GetRunList
    'fixPar','E0 Bkg Norm',...         % free Parameter !!
    'DataType','Real',...              % Real, Twin or Fake
    'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
    'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
    'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
    'chi2','chi2Stat',...              % statistics only
    'NonPoissonScaleFactor',1};
%% build object of MultiRunAnalysis class
obj = MultiRunAnalysis(RunAnaArg{:});

%%
switch obj.DataSet
    case 'Knm2'
        E0str= '_E018573.70eV';
    case 'Knm1'
        E0str= '_E018573.73eV';
end

   
for i = 1:obj.nRuns
    progressbar(i/obj.nRuns);
    filenameT = [strrep(obj.RunData.matFilePath,obj.DataSet,['Twin',obj.DataSet]),'Twin',num2str(obj.RunList(i)),E0str,'.mat'];
    filenameD = [obj.RunData.matFilePath,num2str(obj.RunList(i)),'.mat'];
    %disp(filename);
    
    d   = importdata(filenameD);
    MACE_Bmax_T    = d.MACE_Bmax_T;
    RW_BiasVoltage = d.RW_BiasVoltage;
    save(filenameT,'MACE_Bmax_T','RW_BiasVoltage','-append');
end



