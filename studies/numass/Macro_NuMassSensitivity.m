addpath(genpath('../../../Samak2.0'));

%
Scaling = 100;
[ dr30_b1_off_m dr30_b1_off_p] = NuMassSensitivity('CovMat','OFF','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',1e-3)
[ dr30_b1_on_m dr30_b1_on_p] = NuMassSensitivity('CovMat','OFF','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',1e-3)
[ dr30_b10_off_m dr30_b10_off_p] = NuMassSensitivity('CovMat','OFF','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',10e-3)
[ dr30_b10_on_m dr30_b10_on_p] = NuMassSensitivity('CovMat','OFF','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',10e-3)
[ dr30_b100_off_m dr30_b100_off_p] = NuMassSensitivity('CovMat','OFF','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',100e-3)
[ dr30_b100_on_m dr30_b100_on_p] = NuMassSensitivity('CovMat','OFF','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',100e-3)

fprintf('\n-------------------------------------------\n');
fprintf(' STAT ONLY \n');
fprintf('\n-------------------------------------------\n');

fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 1 mcps - OFF - %g / %g meV \n',dr30_b1_off_m,dr30_b1_off_p);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 1 mcps - ON - %g / %g meV \n',dr30_b1_on_m,dr30_b1_on_p);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 10 mcps - OFF - %g / %g meV \n',dr30_b10_off_m,dr30_b10_off_p);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 10 mcps - ON - %g / %g meV \n',dr30_b10_on_m,dr30_b10_on_p);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 100 mcps - OFF - %g / %g meV \n',dr30_b100_off_m,dr30_b100_off_p);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 100 mcps - ON - %g / %g meV \n',dr30_b100_on_m,dr30_b100_on_p);
fprintf('-------------------------------------------\n');



%
[ dr30_b1_off_m_u dr30_b1_off_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',1e-3)
[ dr30_b1_on_m_u dr30_b1_on_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',1e-3)
[ dr30_b10_off_m_u dr30_b10_off_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',10e-3)
[ dr30_b10_on_m_u dr30_b10_on_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',10e-3)
[ dr30_b100_off_m_u dr30_b100_off_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',100e-3)
[ dr30_b100_on_m_u dr30_b100_on_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',100e-3)

[ flat60_b1_off_m_u flat60_b1_off_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','OFF','BKG_RateAllFPDSec',1e-3)
[ flat60_b1_on_m_u flat60_b1_on_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','ON','BKG_RateAllFPDSec',1e-3)
[ flat60_b10_off_m_u flat60_b10_off_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','OFF','BKG_RateAllFPDSec',10e-3)
[ flat60_b10_on_m_u flat60_b10_on_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','ON','BKG_RateAllFPDSec',10e-3)
[ flat60_b100_off_m_u flat60_b100_off_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','OFF','BKG_RateAllFPDSec',100e-3)
[ flat60_b100_on_m_u flat60_b100_on_p_u] = NuMassSensitivity('CovMat','ON','SystCorr','OFF',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','ON','BKG_RateAllFPDSec',100e-3)

%

[ dr30_b1_off_m_c dr30_b1_off_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',1e-3)
[ dr30_b1_on_m_c dr30_b1_on_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',1e-3)
[ dr30_b10_off_m_c dr30_b10_off_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',10e-3)
[ dr30_b10_on_m_c dr30_b10_on_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',10e-3)
[ dr30_b100_off_m_c dr30_b100_off_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','OFF','BKG_RateAllFPDSec',100e-3)
[ dr30_b100_on_m_c dr30_b100_on_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','DR30','PS_Wein93','ON','BKG_RateAllFPDSec',100e-3)

[ flat60_b1_off_m_c flat60_b1_off_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','OFF','BKG_RateAllFPDSec',1e-3)
[ flat60_b1_on_m_c flat60_b1_on_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','ON','BKG_RateAllFPDSec',1e-3)
[ flat60_b10_off_m_c flat60_b10_off_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','OFF','BKG_RateAllFPDSec',10e-3)
[ flat60_b10_on_m_c flat60_b10_on_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','ON','BKG_RateAllFPDSec',10e-3)
[ flat60_b100_off_m_c flat60_b100_off_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','OFF','BKG_RateAllFPDSec',100e-3)
[ flat60_b100_on_m_c flat60_b100_on_p_c] = NuMassSensitivity('CovMat','ON','SystCorr','ON',...
    'mnuSq_nTrials',Scaling,'pub','ON','display','OFF','Fitter','Minuit',...
    'TD','Flat60','PS_Wein93','ON','BKG_RateAllFPDSec',100e-3)


%
fprintf('\n-------------------------------------------\n');
fprintf(' SYS CORR OFF \n');
fprintf('\n-------------------------------------------\n');

fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 1 mcps - OFF - %g / %g meV \n',dr30_b1_off_m_u,dr30_b1_off_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 1 mcps - ON - %g / %g meV \n',dr30_b1_on_m_u,dr30_b1_on_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 10 mcps - OFF - %g / %g meV \n',dr30_b10_off_m_u,dr30_b10_off_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 10 mcps - ON - %g / %g meV \n',dr30_b10_on_m_u,dr30_b10_on_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 100 mcps - OFF - %g / %g meV \n',dr30_b100_off_m_u,dr30_b100_off_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 100 mcps - ON - %g / %g meV \n',dr30_b100_on_m_u,dr30_b100_on_p_u);
fprintf('-------------------------------------------\n');

fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 1 mcps - OFF - %g / %g meV \n',flat60_b1_off_m_u,flat60_b1_off_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 1 mcps - ON - %g / %g meV \n',flat60_b1_on_m_u,flat60_b1_on_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 10 mcps - OFF - %g / %g meV \n',flat60_b10_off_m_u,flat60_b10_off_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 10 mcps - ON - %g / %g meV \n',flat60_b10_on_m_u,flat60_b10_on_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 100 mcps - OFF - %g / %g meV \n',flat60_b100_off_m_u,flat60_b100_off_p_u);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 100 mcps - ON - %g / %g meV \n',flat60_b100_on_m_u,flat60_b100_on_p_u);
fprintf('-------------------------------------------\n');


%
fprintf('\n-------------------------------------------\n');
fprintf(' SYS CORR ON \n');
fprintf('\n-------------------------------------------\n');

fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 1 mcps - OFF - %g / %g meV \n',dr30_b1_off_m_c,dr30_b1_off_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 1 mcps - ON - %g / %g meV \n',dr30_b1_on_m_c,dr30_b1_on_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 10 mcps - OFF - %g / %g meV \n',dr30_b10_off_m_c,dr30_b10_off_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 10 mcps - ON - %g / %g meV \n',dr30_b10_on_m_c,dr30_b10_on_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 100 mcps - OFF - %g / %g meV \n',dr30_b100_off_m_c,dr30_b100_off_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'DR30 - 100 mcps - ON - %g / %g meV \n',dr30_b100_on_m_c,dr30_b100_on_p_c);
fprintf('-------------------------------------------\n');

fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 1 mcps - OFF - %g / %g meV \n',flat60_b1_off_m_c,flat60_b1_off_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 1 mcps - ON - %g / %g meV \n',flat60_b1_on_m_c,flat60_b1_on_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 10 mcps - OFF - %g / %g meV \n',flat60_b10_off_m_c,flat60_b10_off_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 10 mcps - ON - %g / %g meV \n',flat60_b10_on_m_c,flat60_b10_on_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 100 mcps - OFF - %g / %g meV \n',flat60_b100_off_m_c,flat60_b100_off_p_c);
fprintf('-------------------------------------------\n');
fprintf('\n-------------------------------------------\n');
fprintf(2,'Flat60 - 100 mcps - ON - %g / %g meV \n',flat60_b100_on_m_c,flat60_b100_on_p_c);
fprintf('-------------------------------------------\n');

