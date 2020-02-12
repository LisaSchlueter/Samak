addpath(genpath('../../../Samak2.0'));
runNr = [40257:40260, 40263:40266];
for i=1:numel(runNr)
ft_name = ['spectrum_',num2str(runNr(i))];%,'_outerRingsExcl'];
d = importdata(strcat(ft_name,'.txt'));

qU = d(:,1);
Counts = d(:,2);
CountsErr = d(:,3);
if mod(runNr(i),2)==1 %uneven runs are downwards scanned
    qU = flip(qU);
    Counts = flip(Counts);
    CountsErr = flip(CountsErr);
end
save_name = strcat('../../../samaktest/inputs/FT-Data/',ft_name,'.mat');
save(save_name,'qU','Counts','CountsErr');
end