%% ====    KSN 1 ANALYSIS - SANITY CHECKS    ==== %%

diary 'test_status.txt'
fprintf('\n====================================\n')
fprintf('====    STARTING TEST SERIES    ====\n')
fprintf('====================================\n\n')
diary off



%% ===   STAT+1 TESTS   === %%

diary 'test_status.txt'
fprintf('\n== STARTING STAT+1 ==\n\n')
diary off

%% 1 - RF_EL

SysEffects = struct(...
    'RF_EL','ON',...        % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','OFF',...         % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency

savename   = 'coord_stat+RF_EL.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSTAT + RF_EL DONE (1/9)\n')
diary off

%% 2 - RF_BF

SysEffects = struct(...
    'RF_EL','OFF',...
    'RF_BF','ON',...
    'RF_RX','OFF',...
    'FSD','OFF',...
    'TASR','OFF',...
    'TCoff_RAD','OFF',...
    'TCoff_OTHER','OFF',...
    'DOPoff','OFF',...
    'Stack','OFF',...
    'FPDeff','OFF');

savename   = 'coord_stat+RF_BF.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSTAT + RF_BF DONE (2/9)\n')
diary off


%% 3 - RF_BX

SysEffects = struct(...
    'RF_EL','OFF',...       % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','ON',...        % Column Density, inel cross ection
    'FSD','OFF',...         % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency

savename   = 'coord_stat+RF_BX.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSTAT + RF_BX DONE (3/9)\n')
diary off


%% 4 - FSD
SysEffects = struct(...
    'RF_EL','OFF',...        % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','ON',...          % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency


savename   = 'coord_stat+FSD.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSTAT + FSD DONE (4/9)\n')
diary off

%% 5 - TASR
SysEffects = struct(...
    'RF_EL','OFF',...        % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','OFF',...          % final states
    'TASR','ON',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency


savename   = 'coord_stat+TASR.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSTAT + TASR DONE (5/9)\n')
diary off

%% 6 - TCoff_OTHER
SysEffects = struct(...
    'RF_EL','OFF',...        % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','OFF',...          % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','ON',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency


savename   = 'coord_stat+TCoff_OTHER.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSTAT + TCoff_other DONE (6/9)\n')
diary off

%% 7 - Stack
SysEffects = struct(...
    'RF_EL','OFF',...        % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','OFF',...          % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','ON',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency


savename   = 'coord_stat+Stack.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSTAT + Stack DONE (7/9)\n')
diary off

%% 8 - FPDeff
SysEffects = struct(...
    'RF_EL','OFF',...        % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','OFF',...          % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','ON');        % detector efficiency


savename   = 'coord_stat+FPDeff.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSTAT + FPDeff DONE (8/9)\n')
diary off

%% 9 - BKG
SysEffects = struct(...
    'RF_EL','OFF',...        % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','OFF',...          % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency


savename   = 'coord_stat+BKG.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSTAT + BKG DONE (9/9)\n')
diary off





%% ===   SYST-1 TESTS   === %%

diary 'test_status.txt'
fprintf('\n== STARTING SYST-1 ==\n\n')
diary off

%% 1 - RF_EL

SysEffects = struct(...
    'RF_EL','OFF',...        % Response Function(RF) EnergyLoss
    'RF_BF','ON',...       % RF B-Fields
    'RF_RX','ON',...       % Column Density, inel cross ection
    'FSD','ON',...         % final states
    'TASR','ON',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','ON',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','ON',...       % stacking / HV fluctuations
    'FPDeff','ON');        % detector efficiency

savename   = 'coord_syst-RF_EL.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSYST - RF_EL DONE (1/9)\n')
diary off

%% 2 - RF_BF

SysEffects = struct(...
    'RF_EL','ON',...
    'RF_BF','OFF',...
    'RF_RX','ON',...
    'FSD','ON',...
    'TASR','ON',...
    'TCoff_RAD','OFF',...
    'TCoff_OTHER','ON',...
    'DOPoff','OFF',...
    'Stack','ON',...
    'FPDeff','ON');

savename   = 'coord_syst-RF_BF.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSYST - RF_BF DONE (2/9)\n')
diary off


%% 3 - RF_BX

SysEffects = struct(...
    'RF_EL','ON',...        % Response Function(RF) EnergyLoss
    'RF_BF','ON',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','ON',...         % final states
    'TASR','ON',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','ON',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','ON',...       % stacking / HV fluctuations
    'FPDeff','ON');        % detector efficiency

savename   = 'coord_syst-RF_BX.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSYST - RF_BX DONE (3/9)\n')
diary off


%% 4 - FSD
SysEffects = struct(...
    'RF_EL','ON',...        % Response Function(RF) EnergyLoss
    'RF_BF','ON',...       % RF B-Fields
    'RF_RX','ON',...       % Column Density, inel cross ection
    'FSD','OFF',...          % final states
    'TASR','ON',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','ON',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','ON',...       % stacking / HV fluctuations
    'FPDeff','ON');        % detector efficiency


savename   = 'coord_syst-FSD.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSYST - FSD DONE (4/9)\n')
diary off

%% 5 - TASR
SysEffects = struct(...
    'RF_EL','ON',...        % Response Function(RF) EnergyLoss
    'RF_BF','ON',...       % RF B-Fields
    'RF_RX','ON',...       % Column Density, inel cross ection
    'FSD','ON',...          % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','ON',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','ON',...       % stacking / HV fluctuations
    'FPDeff','ON');        % detector efficiency


savename   = 'coord_syst-TASR.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSYST - TASR DONE (5/9)\n')
diary off

%% 6 - TCoff_OTHER
SysEffects = struct(...
    'RF_EL','ON',...        % Response Function(RF) EnergyLoss
    'RF_BF','ON',...       % RF B-Fields
    'RF_RX','ON',...       % Column Density, inel cross ection
    'FSD','ON',...          % final states
    'TASR','ON',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','ON',...       % stacking / HV fluctuations
    'FPDeff','ON');        % detector efficiency


savename   = 'coord_syst-TCoff_OTHER.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSYST - TCoff_other DONE (6/9)\n')
diary off

%% 7 - Stack
SysEffects = struct(...
    'RF_EL','ON',...        % Response Function(RF) EnergyLoss
    'RF_BF','ON',...       % RF B-Fields
    'RF_RX','ON',...       % Column Density, inel cross ection
    'FSD','ON',...          % final states
    'TASR','ON',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','ON',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','ON');        % detector efficiency


savename   = 'coord_syst-Stack.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSYST - Stack DONE (7/9)\n')
diary off

%% 8 - FPDeff
SysEffects = struct(...
    'RF_EL','ON',...        % Response Function(RF) EnergyLoss
    'RF_BF','ON',...       % RF B-Fields
    'RF_RX','ON',...       % Column Density, inel cross ection
    'FSD','ON',...          % final states
    'TASR','ON',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','ON',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','ON',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency


savename   = 'coord_syst-FPDeff.mat';

sens_func(SysEffects,'ON',savename)

diary 'test_status.txt'
fprintf('\nSYST - FPDeff DONE (8/9)\n')
diary off

%% 9 - BKG
SysEffects = struct(...
    'RF_EL','ON',...        % Response Function(RF) EnergyLoss
    'RF_BF','ON',...       % RF B-Fields
    'RF_RX','ON',...       % Column Density, inel cross ection
    'FSD','ON',...          % final states
    'TASR','ON',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','ON',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','ON',...       % stacking / HV fluctuations
    'FPDeff','ON');        % detector efficiency


savename   = 'coord_syst-BKG.mat';

sens_func(SysEffects,'OFF',savename)

diary 'test_status.txt'
fprintf('\nSYST - BKG DONE (9/9)\n')
diary off

%%

diary 'test_status.txt'
fprintf('\n\n====================================\n')
fprintf('======    TEST SERIES OVER    ======\n')
fprintf('====================================\n\n')
diary off