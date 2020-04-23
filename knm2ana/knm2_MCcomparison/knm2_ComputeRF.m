% calculate response function for MC comparison
SynchrotronFlag = 'ON';
AngTF = 'OFF';
RFConvBinStep = 0.01;
IntMode = 'Conv';
if strcmp(IntMode,'Conv')
    RFStepStr = sprintf('_RFbinStep%.2feV',RFConvBinStep);
elseif strcmp(IntMode,'Integral')
    RFStepStr = '';
end
savename = [getenv('SamakPath'),sprintf('knm2ana/knm2_MCcomparison/results/Response/Knm2_SamakRF%s_%s_Sync%s_ScatTF%s.mat',...
    RFStepStr,IntMode,SynchrotronFlag,AngTF)];

tic;
T = ref_FakeRun_KNM2_RFcomparison('reComputeRF','ON','SynchrotronFlag',SynchrotronFlag,...
    'ISCS','Edep','RFBinStep',RFConvBinStep,'AngularTFFlag',AngTF);
if strcmp(IntMode,'Conv')
    T.RFBinStep = RFConvBinStep;
else
    T.RFBinStep = 0.1;
end
T.InitializeRF('IntMode',IntMode);
toc;

Te           = T.Te;
RF           = T.RF;
RF(RF<1e-50) = 0;
qU           = T.qU;

% save
save(savename,'Te','RF','qU'); % mat file
Write2Txt('filename',strrep(savename,'.mat',''),... % txt file
    'variable',[Te';RF'],...
    'variableName','E Prob',...
    'nCol',2);
