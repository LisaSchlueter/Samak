savedir = [getenv('SamakPath'),'studies/Tim_CrossCheckFSD/'];
tf = [savedir,'FSD_check_2.txt'];
dt = importdata(tf);

% compute samak
A = ref_FakeRun_TimCheck;
A.ComputeTBDDS;
TBDDS = A.TBDDS;

Aoff =  ref_FakeRun_TimCheck('TTFSD','OFF','DTFSD','OFF','HTFSD','OFF');
Aoff.ComputeTBDDS;
TBDDS_off = Aoff.TBDDS;

% check binning
if max(abs(dt(1,:)-A.Te'))>1e-05
    fprintf(2,'not the same binning \n');
    return
end

% save for Tim
Variables = [A.Te,TBDDS_off./simpsons(Aoff.Te,TBDDS_off),TBDDS./simpsons(A.Te,TBDDS)]';
savefileS = sprintf('%sFSD_check_Samak',savedir);
Write2Txt('filename',savefileS,'nCol',3,'variableName','Energy DiffSpecwoFSD DiffSpecwFSD',...
    'variable',Variables)
%% abs
GetFigure;
ps    = plot(A.Te-18575,TBDDS./simpsons(A.Te,TBDDS),'LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
pt    = plot(dt(1,:)-18575,dt(3,:)'./simpsons(dt(1,:),dt(3,:)),'LineWidth',2,'Color',rgb('Orange'));
hold on;
psoff = plot(Aoff.Te-18575,TBDDS_off./simpsons(Aoff.Te,TBDDS_off),':','LineWidth',2,'Color',rgb('DodgerBlue'));
ptoff = plot(dt(1,:)-18575,dt(2,:)'./simpsons(dt(1,:),dt(2,:)),':','LineWidth',2,'Color',rgb('DarkRed'));
PrettyFigureFormat
leg = legend([ps,pt,psoff,ptoff],'FSD Samak','FSD Tim','FSD off Samak','FSD off Tim');
PrettyLegendFormat(leg);
xlabel('Energy - 18575');
ylabel('Diff. spec');
set(gca,'YScale','log');
%print([savedir,'TimCheckAbs.png'],'-dpng','-r300');

%% diff
GetFigure;
p    = plot(A.Te-18575,dt(3,:)'./simpsons(dt(1,:),dt(3,:))-TBDDS./simpsons(A.Te,TBDDS),'LineWidth',2);
hold on;
poff = plot(Aoff.Te-18575,dt(2,:)'./simpsons(dt(1,:),dt(2,:))-TBDDS_off./simpsons(Aoff.Te,TBDDS_off),':','LineWidth',2);

PrettyFigureFormat('FontSize',20);
leg = legend([p,poff],'With FSD','Without FSD','Location','northwest');
PrettyLegendFormat(leg);
xlabel('Energy - 18575');
ylabel('Diff. spec difference (Tim - Samak)');
print([savedir,'TimCheckDiff.png'],'-dpng','-r300');

%% ratio
GetFigure;
Ratio = (dt(3,:)'./simpsons(dt(1,:),dt(3,:)))./(TBDDS./simpsons(A.Te,TBDDS));
p    = plot(A.Te-18575,Ratio,'LineWidth',2);
hold on;
poff = plot(Aoff.Te-18575,(dt(2,:)'./simpsons(dt(1,:),dt(2,:)))./(TBDDS_off./simpsons(Aoff.Te,TBDDS_off)),':','LineWidth',2);

PrettyFigureFormat('FontSize',20);
leg = legend([p,poff],'With FSD','Without FSD','Location','northwest');
PrettyLegendFormat(leg);
xlabel('Energy - 18575');
ylabel('Diff. spec ratio (Tim / Samak)');
ylim([0.89,1.11])
print([savedir,'TimCheckRatio.png'],'-dpng','-r300');