
load('data/R40667SinglePix.mat');

% SP.TableShortStat(SP.TableShortStat(:,13) < 0,:) = [];
% SP.TableShortStat(91,:) = [];

E0_i = 18573.7;


% Converting bias to absolute values

TableShortStat(:,2) = TableShortStat(:,2) + E0_i;
TableShortStat(:,3) = abs(TableShortStat(:,3) + MRA.ModelObj.BKG_RateSec_i);
TableShortStat(:,4) = TableShortStat(:,4) + 1;

TableMedStat(:,2) = TableMedStat(:,2) + E0_i;
TableMedStat(:,3) = abs(TableMedStat(:,3) + MRA.ModelObj.BKG_RateSec_i);
TableMedStat(:,4) = TableMedStat(:,4) + 1;

TableLongStat(:,2) = TableLongStat(:,2) + E0_i;
TableLongStat(:,3) = abs(TableLongStat(:,3) + MRA.ModelObj.BKG_RateSec_i);
TableLongStat(:,4) = TableLongStat(:,4) + 1;

dlmwrite('R40667_Short_Stat.dat',TableShortStat,'precision',10);
dlmwrite('R40667_Med_Stat.dat',TableMedStat,'precision',10);
dlmwrite('R40667_Long_Stat.dat',TableLongStat,'precision',10);



w_SP_s = 1./TableShortStat(:,8).^2;
E_SP_s = wmean(TableShortStat(:,2),w_SP_s);
s_SP_s = 1./sqrt(sum(w_SP_s));

w_SP_m = 1./TableMedStat(:,8).^2;
E_SP_m = wmean(TableMedStat(:,2),w_SP_m);
s_SP_m = 1./sqrt(sum(w_SP_m));

w_SP_l = 1./TableLongStat(:,8).^2;
E_SP_l = wmean(TableLongStat(:,2),w_SP_l);
s_SP_l = 1./sqrt(sum(w_SP_l));