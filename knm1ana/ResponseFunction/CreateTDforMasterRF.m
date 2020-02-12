%% Build Master RF TD
qu = [18373  18483  18488  18493  18498  18503  18508  18513  18518  18523  18528  18533  18535  18537  18539  18541  18543  18545  18547  18549  18551  18553  18555  18557  18559  18561  18562  18563  18564  18565  18566  18567  18569  18571  18573  18578  18583  18593  18603  18623];

qU      = qu';
qUfrac  = 1/numel(qU).*ones(numel(qU),1);
TD      = 'MasterRF';
RunTime = 7200;
save([TD '.mat'],'qU','qUfrac','RunTime','TD');
!mv MasterRF.mat ../../simulation/katrinsetup/TD_DataBank/