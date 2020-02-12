rhodvalues = [126,116,104,135,126,115]-100;
up_unc = [13.17,4.43,3.91,26.24,11.6,7.3];
low_unc = [13.36,5.8,3.13,21.67,12.9,7.2];

fr = 1:3;
sr = 4:6;

hold on
bfr = bar(fr,rhodvalues(fr),'facecolor',rgb('cadetblue'));
errorb(fr,rhodvalues(fr),up_unc(fr),low_unc(fr));
bsr = bar(sr,rhodvalues(sr),'facecolor',rgb('Tomato'));
errorb(sr,rhodvalues(sr),up_unc(sr),low_unc(sr));
hold off
legend([bfr,bsr],{'stat.','stat. + sys.'})
title('CD estimation (100 % = 4.46 \times 10^{17} mol/cm^2)')

ylabel('percentage (relative to run 40668) - 100 %');   
xticks(1:6)
xticklabels(["Short","Med.","Long","Short",...
    "Med.","Long"])
xtickangle(45);

PrettyFigureFormat;

export_fig('plots/rhod_handle.pdf','-pdf');
export_fig('plots/rhod_handle.png','-png');
