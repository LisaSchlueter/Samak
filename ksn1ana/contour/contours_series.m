%% -- KATRIN KSN1 ANALYSIS - SENSITIVITY CONTOUR SERIES -- %%

serie1 = ['FIX','FREE'];
serie2 = ['stat','syst'];
serie3 = [95,90];

for par1 = serie1
    for par2 = serie2
        for par3 = serie3
            % Parameter initialisation
            contour = struct(...
                'CL',par3,...
                'datatype','Twin',...
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