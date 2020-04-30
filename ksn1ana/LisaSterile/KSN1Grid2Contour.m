function [mnu4sq_contour, sin2t4_contour] = KSN1Grid2Contour(mnu4sq,sin2t4,chi2,chi2_ref,CL)
%% calculate contour from msq4 - sin2t4 grid at confidence level cl

% define delta chi2 for given confidence level
DeltaChi2 = GetDeltaChi2(CL/100,2);

nGridStep = size(mnu4sq,1);
sin2t4_contour = zeros(nGridStep,1);
mnu4sq_contour = mnu4sq(:,1);

for i=1:nGridStep
    sin2t4_contour(i) = interp1(chi2(:,i)-chi2_ref,sin2t4(i,:),DeltaChi2,'spline');
end

ExclLogic = abs(sin2t4_contour)>=1;
sin2t4_contour(ExclLogic) = [];
mnu4sq_contour(ExclLogic) = [];
end