% calculate response function for MC comparison
SynchrotronFlag = 'ON';
savename = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/Knm2_SamakRF.mat'];

if strcmp( SynchrotronFlag,'OFF')
    savename = strrep(savename,'.mat','NoSynchrotron.mat');
end

T = ref_FakeRun_KNM2_RFcomparison('reComputeRF','ON','SynchrotronFlag',SynchrotronFlag,'ISCS','Edep');

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
