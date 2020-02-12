
x   = [1 2 3 4];
l   = ['Iterative Fit' 'Multipixel Fit' 'Rest et al' 'M.Slezák'];
e   = [30472.569  30472.564 30472.592 30472.581];
ee  = [0.012 0.014 0.018 0.015];
w   = [1.149 1.143 1.105 1.135];
we  = [0.024 0.041 0.039 0.045];
chi2= [1268 1267.3 1257 1317];
dof = [1278 1198 1198 1198];

% Line Position
figure(1)
subplot(2,1,1)
errorbar(x,(e./e(2)*100-100)*1e4,ee./e(2)*100*1e4,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
xticks([1 2 3 4])
xticklabels({'Iterative Fit' 'Multipixel Fit' 'Rest et al*' 'M.Slezák**'})
ylabel('Line Position (ppm)')
xlim([0.5 4.5]) 
grid on
set(gca,'FontSize',16);

subplot(2,1,2)
errorbar(x,w./w(2)*100-100,we./w(2)*100,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
ylabel('Line Width (%)')
xlim([0.5 4.5]) 
xticks([1 2 3 4])
xticklabels({'Iterative Fit' 'Multipixel Fit' 'Rest et al*' 'M.Slezák**'})
grid on
set(gca,'FontSize',16);

x   = [1 2 3 4 5 6 7 8 ];
l   = ['4057' '4058' '4059' '4060' '4063' '4064' '4065' '4066'];
E0     = [];
E0Err  = [0.012 0.014 0.018 0.015];


% Line Position
fig = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 1400 800]);
subplot(2,1,1)
errorbar(x,(e./e(2)*100-100)*1e4,ee./e(2)*100*1e4,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
xticks([1 2 3 4])
xticklabels({'Iterative Fit' 'Multipixel Fit' 'Rest et al*' 'M.Slezák**'})
ylabel('Line Position (ppm)')
xlim([0.5 4.5]) 
grid on
set(gca,'FontSize',16);

subplot(2,1,2)
errorbar(x,w./w(2)*100-100,we./w(2)*100,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
ylabel('Line Width (%)')
xlim([0.5 4.5]) 
xticks([1 2 3 4])
xticklabels({'Iterative Fit' 'Multipixel Fit' 'Rest et al*' 'M.Slezák**'})
grid on
set(gca,'FontSize',16);