t=importdata('coord_90eV_Twin_stat.mat');
n=importdata('coord_90eV_Twin_stat_95.mat');

figure(1)
ht=loglog(t.sith4_X,t.m4_Y);
hold on
hn=loglog(n.sith4_X,n.m4_Y);
legend([ht hn],'thierry','nathan');
hold off
