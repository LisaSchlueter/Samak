load('./RelicNuBkg/Refs/RefsFromHeizmannPhD.mat');
load('./RelicNuBkg/SensVsMnuSq_KNM1_[0_4].mat');
eta=Sensitivities;
mnuSq=sqrt(ScanPoints);
load('./RelicNuBkg/SensVsMnuSq_TDR_[0_4].mat');
Faes13     = importdata('RelicNuBkg/Refs/Fäs13.txt');
LVV08lower = importdata('RelicNuBkg/Refs/LVV08lower.txt');
LVV08upper = importdata('RelicNuBkg/Refs/LVV08upper.txt');
RW04lower  = importdata('RelicNuBkg/Refs/RW04lower.txt');
RW04upper  = importdata('RelicNuBkg/Refs/RW04upper.txt');
HM05  = [HM05 HM05];
KFM10 = [KFM10 KFM10];
Lob99 = [Lob99 Lob99];
Rob91 = [Rob91 Rob91];
area(linspace(-sqrt(2),2,2),Lob99,1e16,'FaceColor',[0 0.75 0.75]);
hold on;
area(linspace(-sqrt(2),2,2),HM05,1e16,'FaceColor',[0 0.5 0.5]);
area(linspace(-sqrt(2),2,2),Rob91,'FaceColor',[0 0.25 0.25]);
plot(mnuSq,eta,'LineWidth',3,'Color','red');
plot(sqrt(ScanPoints),Sensitivities,'LineWidth',3,'Color','red','LineStyle',':');
plot(linspace(-sqrt(2),2,2),KFM10,'LineWidth',2);
plot(Faes13(:,1),Faes13(:,2),'LineWidth',2,'LineStyle','--');
plot(LVV08lower(:,1),LVV08lower(:,2));
plot(LVV08upper(:,1),LVV08upper(:,2));
x = [LVV08lower(:,1)', fliplr(LVV08upper(:,1)')];
inBetween = [LVV08lower(:,2)', fliplr(LVV08upper(:,2)')];
fill(x, inBetween, [1 1 0]);
plot(RW04lower(:,1),RW04lower(:,2));
plot(RW04upper(:,1),RW04upper(:,2));
x = [RW04lower(:,1)', fliplr(RW04upper(:,1)')];
inBetween = [RW04lower(:,2)', fliplr(RW04upper(:,2)')];
fill(x, inBetween, [0 1 0]);
ax=gca;
ax.YScale='log';
grid on;
text(0.1,1e15,'Robertson et al.','FontSize',14);
text(0.8,5e13,'Hwang et al.','FontSize',14);
text(1.6,1e12,'Lobashev et al.','FontSize',14);
text(1,5e11,'KRN1 Twins','FontSize',14,'FontWeight','bold');
text(0.6,3e10,'KATRIN Sensitivity limit (Bkg rate 130 mcps)','FontSize',14,'FontWeight','bold');
text(0.1,1e9,'Kaboth et al.: KATRIN sensitivity limit (design Bkg)','FontSize',14);
text(1.3,7e6,'Fässler et al.','FontSize',14);
text(1,5e3,'Lazauskas et al.','FontSize',14);
text(1,10,sprintf('Ringwald & Wong \n (extrapolated)'),'FontSize',14);
xlabel('m_{\nu}');
ylabel('\eta (90% C.L.)');
xlim([0 2]);
PrettyFigureFormat;