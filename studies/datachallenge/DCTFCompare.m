clear;
close all;

addpath(genpath('../../../Samak2.0'));

%% Test RF SAMAK and SSC 

% A_S = InitKATRINDC('TD','Flat30','KTFFlag','MACE');
% A_M = InitKATRINDC('TD','Flat30','KTFFlag','MACEMartin');
% RFrange = (18574:0.05:18580)'; % -1 +5 
% 
% samak(:,1) = RFrange;
% martin(:,1) = RFrange;
% 
% samak(:,2) = A_S.KTF(RFrange,18575,A_S.MACE_R_eV);
% martin(:,2) = A_M.KTF(RFrange,18575);

%% Test RF Dom, Florian, Moritz and Woonqook

SSCD = importdata('DC_Data/RSComparison/SSCD.txt');
SSCF = importdata('DC_Data/RSComparison/SSCF.txt');
SSCMW1 = importdata('DC_Data/RSComparison/SSCMW.txt');
SAMAK = importdata('DC_Data/RSComparison/SAMAK.txt');
% 
figure(1)
hold on
plot(SAMAK(:,1),SAMAK(:,2)/trapz(SAMAK(:,1),SAMAK(:,2)));
plot(SSCD(:,1),SSCD(:,2)/trapz(SSCD(:,1),SSCD(:,2)));
plot(SSCF(:,1),SSCF(:,2)/trapz(SSCF(:,1),SSCF(:,2)));
plot(SSCMW1(:,1),SSCMW1(:,2)/trapz(SSCMW1(:,1),SSCMW1(:,2)));

hold off

title(['INTEGRAL spec. normalized',' Comparison '])
xlabel('E - E_0 [eV]');
ylabel('cps per eV');
legend({'SAMAK','SSC Dom','SSC MW','SSC Flo'})

% comparespectra('SAMAKdata',SAMAK,'spec1','SAMAK ','SSCdata',SSCD,'spec2','SSC Dom','compType','INT','errorbarFlag','OFF','normFlag','OFF')
% comparespectra('SAMAKdata',SAMAK,'spec1','SAMAK','SSCdata',SSCMW,'spec2','SSC MW','compType','INT','errorbarFlag','OFF','normFlag','OFF')
% comparespectra('SAMAKdata',SAMAK,'spec1','SAMAK','SSCdata',SSCF,'spec2','SSC Flo','compType','INT','errorbarFlag','OFF','normFlag','OFF')
% comparespectra('SAMAKdata',SSCD,'spec1','SSC Dom','SSCdata',SSCMW,'spec2','SSC MW','compType','INT','errorbarFlag','OFF','normFlag','OFF')
% comparespectra('SAMAKdata',SSCD,'spec1','SSC Dom','SSCdata',SSCF,'spec2','SSC Flo','compType','INT','errorbarFlag','OFF','normFlag','OFF')
% comparespectra('SAMAKdata',SSCF,'spec1','SSC Flo','SSCdata',SSCMW,'spec2','SSC MW','compType','INT','errorbarFlag','OFF','normFlag','OFF')






