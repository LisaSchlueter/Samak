d = importdata('SamakKNM2_DriftInRW123PSR1234_mVperDay.mat');

PlasmaDrift.PseudoRing1 = cell2mat(d.SlopeRW123PSR1234_mV(:,1));
PlasmaDrift.PseudoRing2 = cell2mat(d.SlopeRW123PSR1234_mV(:,2));
PlasmaDrift.PseudoRing3 = cell2mat(d.SlopeRW123PSR1234_mV(:,3));
PlasmaDrift.PseudoRing4 = cell2mat(d.SlopeRW123PSR1234_mV(:,4));
PlasmaDriftErr.PseudoRing1 = cell2mat(d.SlopeErrorRW123PSR1234_mV(:,1));
PlasmaDriftErr.PseudoRing2 = cell2mat(d.SlopeErrorRW123PSR1234_mV(:,2));
PlasmaDriftErr.PseudoRing3 = cell2mat(d.SlopeErrorRW123PSR1234_mV(:,3));
PlasmaDriftErr.PseudoRing4 = cell2mat(d.SlopeErrorRW123PSR1234_mV(:,4));

d2 = importdata('SamakKNM2_ShiftRW123PSR1234_mV.mat');
PlasmaShift.PseudoRing1 = [d2.ShiftRW12PSR1234(1),0,d2.ShiftRW23PSR1234(1)];
PlasmaShift.PseudoRing2 = [d2.ShiftRW12PSR1234(2),0,d2.ShiftRW23PSR1234(2)];
PlasmaShift.PseudoRing3 = [d2.ShiftRW12PSR1234(3),0,d2.ShiftRW23PSR1234(3)];
PlasmaShift.PseudoRing4 = [d2.ShiftRW12PSR1234(4),0,d2.ShiftRW23PSR1234(4)];

PlasmaShiftErr.PseudoRing1 = [d2.ShiftErrorRW12PSR1234(1),0,d2.ShiftErrorRW23PSR1234(1)];
PlasmaShiftErr.PseudoRing2 = [d2.ShiftErrorRW12PSR1234(2),0,d2.ShiftErrorRW23PSR1234(2)];
PlasmaShiftErr.PseudoRing3 = [d2.ShiftErrorRW12PSR1234(3),0,d2.ShiftErrorRW23PSR1234(3)];
PlasmaShiftErr.PseudoRing4 = [d2.ShiftErrorRW12PSR1234(4),0,d2.ShiftErrorRW23PSR1234(4)];

Info = struct('DataSet','Knm2');
Info.PlasmaDrift = 'linear drift over time for each rear wall setting and each pseudo ring';
Info.PlasaShift =  'Shift of mean plasma potential for each rear wall setting with respect to anchor point (RW2, Pseudo-Ring 1) ';
Info.PseudoRings = '4 Pseudo-Rings (RingMerge Full)';
Info.Unit = 'in mV per day';
Info.Method = '300 eV rate monitor analysis';

savedir = [getenv('SamakPath'),'inputs/Plasma/Knm2/'];
MakeDir(savedir);
savename = [savedir,'Knm2_PlasmaTimeEvolution_RMpoint.mat'];

save(savename,'Info',...
    'PlasmaDrift','PlasmaDriftErr',...
    'PlasmaShift','PlasmaShiftErr');

fprintf('save file to %s \n',savename)
