% Plot column density profile
zprofile = (0:1:99)';
profile = [5.16803e+19            8.76746e+19             1.2001e+20            1.49401e+20            1.76413e+20            2.01489e+20            2.24979e+20            2.47154e+20            2.68225e+20            2.88356e+20            3.07673e+20            3.26274e+20            3.44236e+20             3.6162e+20            3.78475e+20             3.9484e+20            4.10751e+20            4.26236e+20            4.41322e+20            4.56034e+20            4.70395e+20            4.84426e+20            4.98148e+20            5.11581e+20            5.24743e+20            5.37651e+20            5.50323e+20            5.62773e+20            5.75016e+20            5.87064e+20            5.98927e+20            6.10617e+20            6.22141e+20  6.335059999999999e+20            6.44719e+20            6.55785e+20            6.66708e+20            6.77492e+20             6.8814e+20  6.986540000000001e+20  7.090380000000001e+20            7.19295e+20            7.29427e+20             7.3944e+20  7.493379999999999e+20  7.591260000000001e+20            7.68812e+20            7.78401e+20            7.87904e+20  7.973299999999999e+20  7.973299999999999e+20            7.87904e+20            7.78401e+20            7.68812e+20  7.591260000000001e+20  7.493379999999999e+20             7.3944e+20            7.29427e+20            7.19295e+20  7.090380000000001e+20  6.986540000000001e+20             6.8814e+20            6.77492e+20            6.66708e+20            6.55785e+20            6.44719e+20  6.335059999999999e+20            6.22141e+20            6.10617e+20            5.98927e+20            5.87064e+20            5.75016e+20            5.62773e+20            5.50323e+20            5.37651e+20            5.24743e+20            5.11581e+20            4.98148e+20            4.84426e+20            4.70395e+20            4.56034e+20            4.41322e+20            4.26236e+20            4.10751e+20             3.9484e+20            3.78475e+20             3.6162e+20            3.44236e+20            3.26274e+20            3.07673e+20            2.88356e+20            2.68225e+20            2.47154e+20            2.24979e+20            2.01489e+20            1.76413e+20            1.49401e+20             1.2001e+20            8.76746e+19            5.16803e+19]';
z = (0:0.1:99)';
CD = interp1(zprofile,profile./(5e21),z,'spline');

f3 = figure('Renderer','opengl');
set(f3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7 ,0.7]);
plot(z./100,CD,'LineWidth',5,'Color',rgb('CadetBlue'));
PrettyFigureFormat;
set(gca,'FontSize',24);
ylim([min(CD),max(CD)*1.02]);
ylabel('rel. column density');
xlabel('rel. position z in WGTS');

fig_str = 'ColumnDensityProfile';
publish_figurePDF(gcf,['../ft_CovarianceMatrices/plots/pdf/',fig_str,'.pdf']);
print(gcf,['../ft_CovarianceMatrices/plots/png/',fig_str,'.png'],'-dpng','-r400');
