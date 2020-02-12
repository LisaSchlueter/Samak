%
% Read / Load Magnetic Field as a function of z
% Interpolate B field as a function of z
%
function b = bfield(z)
global LookUpz;
global LookUpB;

b = interp1(LookUpz,LookUpB,z);
end
