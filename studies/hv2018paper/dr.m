
x   = [1 2 3 4 5 6 7 8] ;
l   = ['40257' '40258' '40259' '40260' '40263' '40264' '40265' '40266'];
e   = [0.295 0.319 0.353 0.268 0.350 0.308 0.311 0.316]*1000;
ee  = [0.023 0.024 0.021 0.024 0.023 0.023 0.023 0.023]*1000;


figure(1)
errorbar(x,e,ee,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
xticks([1 2 3 4 5 6 7 8])
xticklabels({'40257' '40258' '40259' '40260' '40263' '40264' '40265' '40266'})
ylabel('Background Rate (mcps)')
xlim([0.5 8.5]) 
grid on
set(gca,'FontSize',16);
