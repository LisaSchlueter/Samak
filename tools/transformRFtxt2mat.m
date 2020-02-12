clear;
addpath(genpath('../../../Samak2.0'));

for qU = [18545 18560 18575]
    ssc  = importdata(['../../inputs/TF/RF_ssc_nScat6_DetailedOff_Xsection3p46_Bs3p6T_Ba5p83051G_BMax6T_ColumnDensity_5E21_at',num2str(qU),'.txt']);
    save(['../../inputs/TF/RF_ssc_nScat6_DetailedOff_Xsection3p46_Bs3p6T_Ba5p83051G_BMax6T_ColumnDensity_5E21_at',num2str(qU),'.mat'],'ssc','-v7.3','-nocompression')
end

% for qU = 1:1
% bckstruct  = importdata(['../../inputs/BCK/BackgroundRate_35023-35110','','.txt']);
% bck = bckstruct.data;
% save(['../../inputs/BCK/BackgroundRate_35023-35110','','.mat'],'bck','-v7.3','-nocompression')
% end
%  
% bkc_x = load(['../../inputs/BCK/BackgroundRate_35023-35110','','.mat']);
% qU = 18545;
% ssc_x = load(['../../inputs/TF/RFSSC/RF_ssc_nScat6_DetailedOff_Xsection3p46_Bs3p6T_Ba9G_BMax6T_qU',num2str(qU),'.mat'],'ssc');
