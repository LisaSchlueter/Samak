H0 = importdata('RelicNuBkg/eta0.txt');
H0(:,1)=H0(:,1)+18570;
IntHeiz = simpsons(H0(:,1),H0(:,2));
T0 = ref_RelicNuBkg_DesignReport('ToggleRelic','OFF','mnuSq_i',1,'TimeSec',1);
T0.ComputeTBDDS;
T0.TBDDS = T0.TBDDS.*0.008077/5e17;
IntTDR = simpsons(T0.Te,T0.TBDDS);
Window=(T0.Te>=H0(1,1))&(T0.Te<=H0(end,1));
TBDDS = T0.TBDDS(Window);
Te = T0.Te(Window);
IntTail = simpsons(Te,TBDDS);
Frac=IntTail/IntTDR*T0.CumFrac;
h0 = semilogy(H0(:,1)-T0.Q,H0(:,2),'LineWidth',2,'Color','Red','LineStyle','--');
hold on;
t0 = semilogy((Te-T0.Q),TBDDS,'LineWidth',2,'Color','Blue','LineStyle','-.');
grid on;
ylim([1e-23,1e-18]);
xlim([-5,2]);
PrettyFigureFormat;
hold off;
Int = IntHeiz/Frac
int = IntTDR/T0.CumFrac