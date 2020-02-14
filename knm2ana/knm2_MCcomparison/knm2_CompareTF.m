

% Samak TF without Synchrotron
A = ref_FakeRun_KNM2_RFcomparison('SynchrotronFlag','OFF');
qu = 18545;
te = (0:0.01:3)+qu;
MaceTF = A.ComputeMaceTF(te,qu);
plot(te-qu,MaceTF);
savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
savename = [savedir,'Samak_TF'];
Write2Txt('filename',savename,'nCol',2,'variable',[te;MaceTF]);