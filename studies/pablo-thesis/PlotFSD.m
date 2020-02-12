figfsd = figure(10);

obj = ref_RunAnalysis('FAKEFTALLex2','','');
type = 'log';

pos = get(figfsd,'position');
set(figfsd,'position',[pos(1:2)/3.5 pos(3)*1.5,pos(4)*1])
hold on
h1 = stairs(-obj.HTexE_G+1.7,obj.HTexP_G*100,'Color','Black','LineWidth',2);
hold on
h2 = stairs(-obj.HTexE_E+1.7,obj.HTexP_E*100,'Color','Red','LineWidth',2);
hold off
set(gca, 'XScale', 'lin');
switch type
    case 'lin'
        set(gca, 'YScale', 'lin');
    case 'log'
        set(gca, 'YScale', 'log')  ;
        %set(gca, 'XScale', 'log');
end
grid on
%axis([0.1 250 1e-3 50]);
ylim([1e-3 50]);
xlim([-250 10]);
xlabel('Binding Energy (eV)','FontSize',14);
ylabel('Probability (%)','FontSize',14);
fsHTitle=sprintf('HT Ground and Excited States - %s',obj.HTFSD);
title(fsHTitle,'FontSize',14);
set(gca,'FontSize',12);
strG = sprintf('Ground States: (P=%.1f%%)',obj.HTNormGS*100);
strE = sprintf('Excited States: (P=%.1f%%)',obj.HTNormES*100);
a = legend([h1 h2],strG,strE,'Location','NorthWest'); %legend(a,'boxoff','FontSize',12);
PrettyFigureFormat;

export_fig('plots/HTFSD.pdf','-pdf')
