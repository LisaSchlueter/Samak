clear all
range=60;
c=0;
for ba=3:1:9
c=c+1;
[qU, qUfrac, TD]=MTDcreator('MACE_Ba_T',ba*1e-4,'AnchorBkg6G',0.350,'Range',range);
TabqU(c,:)      = qU;
TabqUfrac(c,:)  = qUfrac;
end
figure(1)
stairs(TabqUfrac')
xlabel('qU (V)');
ylabel('Fraction of Time');
set(gca,'YScale','lin');
PrettyFigureFormat;
xlwrite(sprintf('MTD_IsoStatPlus350mcps3to9G_%geV.xls',range),[TabqU;TabqUfrac]);

