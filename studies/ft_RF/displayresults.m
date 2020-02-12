
x   = [1 2 3 4 5 6 7 8 9 10 11];
l   = ['5%/1% C' '5%/1% U' '-' '50%/1% C' '50%/1% U' '-' '5%/0.1% C' '5%/0.1% U' '-' '5%/0.1% C' '5%/0.1% U'];
bias   = [-0.02 0.03 0 -0.03 0.07 0 -0.27 -0.24 0 1 0.4];
error  = [0.3 0.8 0 0.4 1.01 0 2.01 3.01 0 11 10];


% Line Position
figure(1)
errorbar(bias,error,'s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',2)
xticks([1 2 3 4 5 6 7 8 9 10 11])
xticklabels({'5%/1% C' '5%/1% U' '-' '50%/1% C' '50%/1% U' '-' '5%/0.1% C' '5%/0.1% U' '-' '5%/0.1% C' '5%/0.1% U'})
ylabel('E0-18575 (eV)')
xlim([0.5 12.5]) 
grid on
set(gca,'FontSize',16);

