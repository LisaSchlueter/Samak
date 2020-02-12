function [par, err, chi2min, errmat] = pminuit(funcstr,Args)

[par, err, chi2min, errmat] = fminuit(funcstr,Args{:});

end