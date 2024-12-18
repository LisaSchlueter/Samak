function [MasterTe MasterRF] = readMasterRF()
% Read Reference Knm1 Response Function
% te: kinetic energy
% masterRF: matrix providing rf value
%           column = pixel
%           line   = qu
% 
rfdir = '../../../responsefunction-master/';
% Retarding Potentials
qu = [18373  18483  18488  18493  18498  18503  18508  18513  18518  18523  18528  18533  18535  18537  18539  18541  18543  18545  18547  18549  18551  18553  18555  18557  18559  18561  18562  18563  18564  18565  18566  18567  18569  18571  18573  18578  18583  18593  18603  18623];

ElectricPotentialCorrection =  18600 + [-18598.0384     -18598.0396       -18598.04     -18598.0385     -18598.0444     -18598.0459      -18598.047     -18598.0475     -18598.0478     -18598.0483      -18598.049     -18598.0492     -18598.0481      -18598.046      -18598.044     -18598.0435     -18598.0559     -18598.0583      -18598.059     -18598.0589     -18598.0591     -18598.0604     -18598.0619     -18598.0619      -18598.059      -18598.055     -18598.0526     -18598.0531     -18598.0643     -18598.0686     -18598.0703     -18598.0699     -18598.0695     -18598.0702     -18598.0728     -18598.0748     -18598.0732     -18598.0679     -18598.0629     -18598.0617     -18598.0772     -18598.0814     -18598.0812     -18598.0801     -18598.0798     -18598.0819     -18598.0858      -18598.087     -18598.0825      -18598.075     -18598.0709      -18598.072     -18598.0841      -18598.091     -18598.0929     -18598.0908     -18598.0897     -18598.0903     -18598.0948      -18598.099     -18598.0977     -18598.0892     -18598.0815     -18598.0799     -18598.0981     -18598.1039     -18598.1023        -18598.1     -18598.0991      -18598.102     -18598.1087     -18598.1114     -18598.1053     -18598.0942     -18598.0889     -18598.0903      -18598.103     -18598.1128     -18598.1143     -18598.1105     -18598.1089     -18598.1089     -18598.1156     -18598.1223     -18598.1214     -18598.1096     -18598.0995     -18598.0976     -18598.1183     -18598.1256      -18598.122     -18598.1189     -18598.1173     -18598.1207     -18598.1303     -18598.1348     -18598.1276     -18598.1126     -18598.1066     -18598.1078      -18598.121      -18598.134     -18598.1348     -18598.1288     -18598.1273     -18598.1261     -18598.1354     -18598.1446     -18598.1445     -18598.1295     -18598.1168     -18598.1148     -18598.1376     -18598.1466     -18598.1406     -18598.1366     -18598.1344     -18598.1381      -18598.151     -18598.1574     -18598.1492     -18598.1299     -18598.1236     -18598.1244      -18598.138     -18598.1545     -18598.1543     -18598.1459     -18598.1448     -18598.1419      -18598.154     -18598.1658      -18598.167     -18598.1483     -18598.1333     -18598.1314     -18598.1561     -18598.1668      -18598.158     -18598.1535     -18598.1504     -18598.1543     -18598.1707      -18598.179     -18598.1704     -18598.1462     -18598.1403       -18598.14];

% Reading Master RF data
% Each file has a name knm1_response_[integer_qU_value].000000.dat
% The file contains the values of response functions for all 148 pixels.
% The first column in the file is the kinetic energy of the electron in eV.
% The second column is RF for given kinetic energy for pixel 0, third - RF
% for pixel 1 and so on.
MasterTe       = zeros(3000,148);
MasterRF       = zeros(numel(qu),3000,148);
counterbis=0;
for i=1:numel(qu)
    masterRFfile = [rfdir 'knm1_response_' num2str(qu(i)) '.000000.dat'];
    rftmp = importdata(masterRFfile);
    disp(masterRFfile);
    % read Kinetic Energy from first file
    if i == 1
        for p=1:148
        MasterTe(:,p) = rftmp(:,1) + ElectricPotentialCorrection(p);
        end
    end
    % read pixelwise response functions
    MasterRF(i,:,:) =  rftmp(:,2:149);
    if (i==1 || i==10 || i==18 || i==40)
    counterbis=counterbis+1;
    MasterRF(counterbis,:,:) =  rftmp(:,2:149);
    end
end
save('MasterRF.mat','MasterTe','MasterRF');
save('MasterRFqU18545.mat','MasterTe','MasterRF');
end
