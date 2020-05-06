%% -- KATRIN KSN1 ANALYSIS - SENSITIVITY CONTOUR -- %%

% Contour parameters
contour = struct(...
    'CL',95,...
    'datatype','Real',...
    'uncertainty','syst',...
    'NPfactor',1,...
    'scan_step',10,...
    'eVrange',90,...
    'activeFlag','OFF');

% Scan
sens_func(contour)