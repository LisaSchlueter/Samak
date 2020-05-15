%% -- KATRIN KSN1 ANALYSIS - SENSITIVITY CONTOUR SERIES -- %%

serie1 = {'OFF'};
serie2 = {'stat','syst'};
serie3 = [40];

for par1 = serie1
    for par2 = serie2
        for par3 = serie3
            % Parameter initialisation
            contour = struct(...
                'CL',par3,...
                'datatype','Real',...
                'uncertainty',par2,...
                'NPfactor',1,...
                'scan_step',10,...
                'eVrange',90,...
                'activeFlag',par1);

            % Scan
            sens_func(contour)
        end
    end
end