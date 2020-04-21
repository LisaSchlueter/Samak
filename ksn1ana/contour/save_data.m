%% Datasave

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,'coord_0to4_d400_40eV.mat'];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');