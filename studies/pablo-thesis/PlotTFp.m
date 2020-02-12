
% A = ref_corrections('FPD_SEGMENTATION','MULTIPIXEL','FPD_Pixel',1:148,...
%     'qUmin',18575-30,'qUmax',18580,'Mode','DataTBD','TD','StackCD100allmpix',...
%     'nTeBinningFactor',50);

B = ref_corrections('FPD_SEGMENTATION','OFF','FPD_Pixel',1:148,...
    'qUmin',18575-30,'qUmax',18580,'Mode','Sim',... %'TD','StackCD100allmpix',...
    'nTeBinningFactor',50,'recomputeRF','OFF');

%load('ResponseFunction/samakRF_5e+17cm2_NIS11_Bm6T_Bs3.6T_Ba3G_multipix.mat');

%RF148 = squeeze(RFip(:,30,:));
% RFp = squeeze(A.RF(:,1,:));
% RF148 = RFp./B.RF(:,1);
plot(B.Te-18545-4,B.RF(:,5));

%axis([0 2 0.95 1.02]);

title('KATRIN Response Function in SAMAK');
ylabel('Transmission Probability');
xlabel('E - qU (eV)');

PrettyFigureFormat;

export_fig('plots/SAMAKRF.png','-png');

