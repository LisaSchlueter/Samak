%
%            FIT INTEGRAL SPECTRUM
%              Fit a 83mKr Lines
%
%                 All Pixels
% 
%          Th. Lasserre - CEA Saclay
%                August 2017
%

close all
ParAllPixels = zeros(148,12);

%TD = 'KrK32';
TD = 'KrL3_32';
%TD = 'KrL3_32_HS';

for p=1:1:10
%for p=1:1:1
    fprintf(2,'**************************\n')
    fprintf(2,'Fitting Pixel %0.2f \n',p)
    [par err chi2min ndof] = KrLineFit('TD',TD,...
        'mypub','OFF','fign',1,...
        'FPD_Segmentation','PIXEL','Pixel',p,...
        'DopplerEffectFlag','ON',...
        'Mode','Data','CPS','OFF',...
        'HVRipples','ON','HVRipplesP2PV',0.52);
    ParAllPixels(p,:) = [(par) err chi2min ndof];
    fprintf(2,'**************************\n')
    disp(ParAllPixels(p,:));
    fprintf(2,'**************************\n')
end

switch TD
    case 'KrK32'
        save('./results/Data_KrK32AllPixelsDopplerONRipplesON.mat','ParAllPixels');
    case 'KrL3_32'
        save('./results/Data_KrL332AllPixelsDopplerONRipplesON.mat','ParAllPixels');
case 'KrL3_32_HS'
        save('./results/Data_KrL332HSAllPixelsDopplerONRipplesON_v3.mat','ParAllPixels');
end
