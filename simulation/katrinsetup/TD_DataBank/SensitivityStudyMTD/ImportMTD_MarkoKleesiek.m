
range = 30;
Ba = [3,6,9];

for i=1:numel(Ba)
file_path = sprintf('../../simulation/katrinsetup/TD_DataBank/SensitivityStudyMTD/MarcoKleesiek/');
file_name = sprintf('mtd_%.0f_350mcps+%.0fG.txt',range,Ba(i));

d = importdata([file_path,file_name]);
qU = d(:,1)+18575.0;
qUfrac = d(:,2)./sum(d(:,2));
TD = sprintf('MK_Sensitivity_BKG350mcps_%.0feV_Ba%.0fG',range,Ba(i));

save(['../../simulation/katrinsetup/TD_DataBank/',TD,'.mat'],'qU','qUfrac','TD',...
    '-v7.3','-nocompression')
end