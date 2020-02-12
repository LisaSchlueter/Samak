% tic; M = MultiRunAnalysis('RunList','KNM1_175mvRW'); toc;
% 
% M.ComputeCM('StackCM','OFF','SysEffect',struct('RF_EL','ON'),'nTrials',100);
% M.FitCM_Obj.SysEffect.RF_BF = 'ON';
% M.FitCM_Obj.SysEffect.RF_RX = 'ON';
%%
% tic;
% M.FitCM_Obj.RecomputeFlag = 'ON';
% M.FitCM_Obj.ComputeCM_RF;
% toc;
%% 
% A = ref_katrin; A.ComputeTBDDS; A.ComputeTBDIS;
% TBDIS1 = A.TBDIS;
% 
% A.WGTS_MolFrac_TT = 0.7;
%  A.AdjustMolFrac; A.ComputeTBDDS;A.ComputeTBDIS;
% TBDIS2 = A.TBDIS;
% 
% B = ref_katrin('WGTS_MolFrac_TT',0.7); B.ComputeTBDDS; B.ComputeTBDIS;
% TBDIS3 = B.TBDIS;
% 
%  plot(A.qU,TBDIS2-TBDIS3);


                    
                    
                    
                    
                    
                    
                    
                    
         