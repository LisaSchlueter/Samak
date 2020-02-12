if ~exist('SPR_R','var') || ~exist('SPR_T','var')
    [SPR_R , SPR_T] = Knm1RealTwin_Create();
end

[RT_parqU, RT_errqU, RT_chi2qU, RT_dofqU] = Knm1RealTwin_qUScan(SPR_R,SPR_T);