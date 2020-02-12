
range = 60;
Ba = 3:1:12;

for i=1:numel(Ba)
file_path = sprintf('../../simulation/katrinsetup/TD_DataBank/SensitivityStudyMTD/MTD%.0feV/',range);
file_name = sprintf('mtd_%.0fE-4.txt',Ba(i));

d = importdata([file_path,file_name]);
qU = d(:,1)+18575.0;
qUfrac = d(:,2)./sum(d(:,2));
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,Ba(i));

save(['../../simulation/katrinsetup/TD_DataBank/',TD,'.mat'],'qU','qUfrac','TD',...
    '-v7.3','-nocompression')
end