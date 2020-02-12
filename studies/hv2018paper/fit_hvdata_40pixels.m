% 
% Fit HV Paper Data 
% L3-32 Line
% 40 pixels (3 inner Rings)
% 
% Th. Lasserre - CEA Saclay
% January 2018
%

close all
format long

Par40Pixels = zeros(40,11);

for p=1:1:40
    fprintf(2,'--------------------------------------------------------------\n')
    fprintf(2,'Fitting Pixel %0.2f \n',p)
    [parA, err, chi2min, ndof] =  fit_hvdata_pixel('Pixel',p,'display','OFF');
    Par40Pixels(p,:) = [p parA err chi2min ndof];
    fprintf(2,'--------------------------------------------------------------\n')
    disp(Par40Pixels(p,:));
    fprintf(2,'--------------------------------------------------------------\n')
end
save('./results/hvdata_l332_40pixels_v3.mat','Par40Pixels');

% Plot the results