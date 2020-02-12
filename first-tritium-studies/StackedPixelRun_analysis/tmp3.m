g   = load('StackCD100allex2_Gradient.mat');
g2  = load('StackCD100allex2_Gradient2.mat');
ng  = load('StackCD100allex2.mat'); 

figure(55)
plot(g.qU,(g.TBDIS - ng.TBDIS)./sqrt(ng.TBDIS),'LineWidth',2);
