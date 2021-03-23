
Iso = 'DT';
path = [getenv('SamakPath'),'inputs/FSD/'];

d = importdata([path,'/FSD_KNM2_',Iso,'_0p1eV.txt']);

Etmp = d(:,1);
E = Etmp(Etmp<45);

Probtmp =  d(:,2);
Prob = Probtmp(Etmp<45);

Write2Txt('filename',[path,'FSD_KNM2_',Iso,'_0p1eV_cut40eV'],'nCol',2,'variable',[E';Prob']);