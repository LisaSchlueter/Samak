
Ba = 1e-04.*(3:12);
Bs = 6; 
WGTS_B_T = 3.6;
BKG = GetBackground('MACE_Ba_T',Ba,'WGTS_B_T',WGTS_B_T);
BKG70 = GetBackground('MACE_Ba_T',Ba,'WGTS_B_T',WGTS_B_T*0.7);

R             = 18575*Ba/Bs;     %Resolution 
R70           = 18575*Ba/(Bs*0.7);

plot(Ba*1e4,R,'--o','LineWidth',3,'Color',rgb('IndianRed'));
hold on;
plot(Ba*1e4,R70,'--o','LineWidth',3,'Color',rgb('CadetBlue'));
xlabel('B_a (10^{-4} T)');
ylabel('Resolution \DeltaE (eV)');
set(gca,'XScale','log');
PrettyFigureFormat;
set(gca,'FontSize',16)
legend('100% magnets', '70% magnets');
grid on;