function [mnu4sq_contour, sin2t4_contour] = GetSterileContour(mnu4sq,sin2t4,chi2,chi2_ref,CL)

%% get contour at X sigma

% Confidence level
switch CL
    case 90
        DeltaChi2 = 4.61;
    case 95
        DeltaChi2 = 5.99;
    case 99
        DeltaChi2 = 9.21;
end

nGridStep = size(mnu4sq,1);

sin2t4_contour = zeros(nGridStep,1);
mnu4sq_contour = mnu4sq(:,1);

%chi2_ref =min(min(chi2));
for i=1:nGridStep
    
    %if any((chi2(:,i)-chi2_ref)<0)
       
   % else
    sin2t4_contour(i) = interp1(chi2(:,i)-chi2_ref,sin2t4(1,:),DeltaChi2,'spline');
    %end
end

end