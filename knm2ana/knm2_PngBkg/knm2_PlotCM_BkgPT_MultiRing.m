

dir = [getenv('SamakPath'),'inputs/CovMat/BkgPT/CM/'];
fname = [dir,'BkgPT_KNM2_Prompt_B220mcps_Slope3.0muCpsPerS_SlopeErr3.0muCpsPerS_50000Trials_RingFull.mat'];
d = importdata(fname);
%%
d.obj.SysEffect.BkgPT = 'ON';
d.obj.SysEffect.RF_BF = 'OFF';d.obj.SysEffect.RF_EL = 'OFF'; d.obj.SysEffect.RF_RX = 'OFF'; d.obj.SysEffect.FSD = 'OFF';d.obj.SysEffect.TASR = 'OFF'; d.obj.SysEffect.TCoff_OTHER = 'OFF'; d.obj.SysEffect.Stack = 'OFF'; d.obj.SysEffect.FPDeff = 'OFF';d.obj.SysEffect.LongPlasma = 'OFF';
d.obj.PlotCM('CovMatInput',d.obj.CovMat,'Mode','CM','qUWindowIndexMax',-135,'savename',sprintf('_Err_%.1fmuCpsPerS_CMtot',d.obj.BKG_PtSlopeErr*1e6),'saveplot','ON');
d.obj.PlotCM('CovMatInput',d.obj.CovMatFracShape,'Mode','Shape','qUWindowIndexMax',-135,'savename',sprintf('_Err_%.1fmuCpsPerS_CMFracShape',d.obj.BKG_PtSlopeErr*1e6),'saveplot','ON');
d.obj.PlotCorr('qUWindowIndexMax',-135,'savename',sprintf('_Err_%.1fmuCpsPerS',d.obj.BKG_PtSlopeErr*1e6),'qUWindowIndexMin',40,'saveplot','ON')
