addpath(genpath('../../../Samak2.0'));
%Init
MyObj = InitKrKATRIN_krl332();
MyObj.HVRipplesP2PV = 0.44;
MyObj.FPD_Pixel = 1;
W_bestfit = 1.557;%1.559%1.512;  %Width for Best Fit with Covariance Matrix [eV] Data
WErr_bestfit = 0.3;%0.16  %0.314; %1 Sigma on Width (Best Fit) [eV] Data
nsteps = 35;  %total number of steps = 2*nsteps+1
mychi2min = zeros(2*nsteps,1);
mychi2min_NoCovMat = zeros(2*nsteps,1);
myWidth = zeros(2*nsteps,1);

%Fit Width in KrLineFitData_CovMat before Fitting!
tic
progressbar('Scanning Width');
for i=1:(2*nsteps+1)
progressbar('(i+nsteps)/(2*nsteps)');
myWidth(i) = W_bestfit + (i-nsteps-1)*(WErr_bestfit/nsteps); %was 3*WErr_...
MyObj.L3_32_W_i= myWidth(i);
[par, err, chi2min, ndof, errmat] = KrL3LineFitData_CovMat('CovMat','ON',...
                   'CovMatCompute',1, 'KrObject', MyObj, 'plots', 'OFF'); 
mychi2min(i) = chi2min;  
[par2, err2, fit_xi2min] = KrL3LineFitData_NoCovMat('KrObject', MyObj, 'plots', 'OFF'); 
mychi2min_NoCovMat(i)=fit_xi2min;       
end       
toc

%% Plot
orange=1/255*[255,140,0];  
%set(gcf,'Color','w') %background white
figure(123);
plot(myWidth*1e3, mychi2min, 'Color', orange, 'LineWidth', 2.5); hold on;
plot(myWidth*1e3, mychi2min_NoCovMat, 'Color', 'blue', 'LineWidth', 2.5);
plot([W_bestfit*1e3 W_bestfit*1e3], [0 9], 'k--', 'LineWidth', 2) %Ref Line
plot([1348.4 1348.4], [0 9], 'k:'); plot([1763.8 1763.8], [0 9], 'k:');
plot([1457.5 1457.5], [0 9], 'k:'); plot([1659.86 1659.86], [0 9], 'k:');
plot([1348.4 1763.8], [9 9] ,'k:');
WerrCovMat = sprintf('W = (1557 - 208.6 + 206.8) eV');
WerrNoCovMat = sprintf('W = (1557 - 99.5 + 102.9) eV');
dim = [.375 .75 .16 .16]; dim2 = [.378 .7 .16 .16];
delete(findall(gcf,'type','annotation'));
a1 = annotation('textbox', dim, 'String', WerrCovMat,'FitBoxToText','on');
a2 = annotation('textbox', dim2, 'String', WerrNoCovMat,'FitBoxToText','on');
a1.FontSize= 15; a1.Color = orange; a1.LineStyle = 'none';
a2.FontSize= 15; a2.Color = 'blue'; a2.LineStyle = 'none';
hold off;
%plot([(W_bestfit+WErr_bestfit)*1e3 (W_bestfit+WErr_bestfit)*1e3], [0 1300], 'k--', 'LineWidth',1); %Ref Line
yticks([0:1:50]); xlabel('Width [eV]', 'FontSize', 20); ylabel('$\mathbf{\chi^2_{min}}$ (27 ndof)', 'Interpreter', 'latex', 'FontSize', 22);
xticks([1200:50:1900]);
xlim([(W_bestfit-WErr_bestfit)*1e3 (W_bestfit+WErr_bestfit)*1e3]);
%ylim([min(mychi2min) max(mychi2min)]);
ylim([0 12]);
%grid on;
%title('Fit to Simulation (Pixel 1) uniform CovMat 2V');
reflineleg = sprintf('W_{bf}=%.0f eV', W_bestfit*1e3);
l1 = legend('CovMat', 'NoCovMat', reflineleg, 'Location', 'southeast');
l1.FontSize = 16;
%export_fig './plots/Comparison/fig/ScanWidthSim_nofluct_Pixel1_whole_uniform2V_lines.fig'
%export_fig './plots/Comparison/ScanWidthSim_nofluct_Pixel1_whole_uniform2V_linestalk.pdf'

figure(124);
% W_bestfit = 1.547;  %myWidth(22);
% upperWidth = myWidth(22:end);
% lowerWidth = myWidth(1:22)+2*(W_bestfit-myWidth(1:22)); 
% lowerchi2min = mychi2min(1:22);
% upperchi2min = mychi2min(22:end);
% lowerchi2min_NoCovMat = mychi2min_NoCovMat(1:22);
% upperchi2min_NoCovMat = mychi2min_NoCovMat(22:end);
lowerWidth = myWidth(1:nsteps)+2*(W_bestfit-myWidth(1:nsteps)); %flipped
upperWidth = myWidth(nsteps+1:end);
lowerchi2min = mychi2min(1:nsteps);
upperchi2min = mychi2min(nsteps+1:end);
lowerchi2min_NoCovMat = mychi2min_NoCovMat(1:nsteps);
upperchi2min_NoCovMat = mychi2min_NoCovMat(nsteps+1:end);


plot(lowerWidth*1e3, lowerchi2min, 'Color', orange, 'Linestyle' , '--'); hold on;
plot(upperWidth*1e3, upperchi2min, 'Color', orange, 'LineWidth', 1);
plot(lowerWidth*1e3, lowerchi2min_NoCovMat, 'Color', 'blue','Linestyle', '--'); 
plot(upperWidth*1e3, upperchi2min_NoCovMat, 'Color', 'blue', 'LineWidth', 1); 
plot([1557+99.5 1557+99.5], [0 9], 'k:'); plot([1763.8 1763.8], [0 9], 'k:');
plot([1557+208.6 1557+208.6], [0 9], 'k:'); plot([1659.86 1659.86], [0 9], 'k:');
plot([1348.4 1763.8], [9 9] ,'k:');
%plot([(W_bestfit+WErr_bestfit)*1e3 (W_bestfit+WErr_bestfit)*1e3], [0 1300], 'k--', 'LineWidth',1); %Ref Line
%plot([(W_bestfit+2*WErr_bestfit)*1e3 (W_bestfit+2*WErr_bestfit)*1e3], [20 35], 'k--', 'LineWidth',1); %Ref Line
yticks([0:1:27]);
grid on;
xlabel('Width [eV]', 'FontSize', 20); ylabel('$\mathbf{\chi^2_{min}}$ (27 ndof)', 'Interpreter', 'latex', 'FontSize', 22);
CovMatleg = sprintf('CovMat \n W < W_{bf}'); NoCovMatleg = sprintf('No CovMat \n W < W_{bf}');
l2 = legend(CovMatleg ,'W > W_{bf}',NoCovMatleg, 'W > W_{bf}', 'Location', 'northwest');
l2.FontSize = 14;
%xlim([(W_bestfit)*1e3 (W_bestfit+WErr_bestfit)*1e3]);
%ylim([min(mychi2min) max(mychi2min)]);
xlim([1550 1800]);
ylim([0 12]);
%title('Fit to Simulation (Pixel 1) uniform CovMat 2V');
hold off;
%export_fig './plots/Comparison/fig/ScanWidthSim_nofluct_Pixel1_uniform2V_lines.fig'
%export_fig './plots/Comparison/ScanWidthSim_nofluct_Pixel1_uniform2V_linestalk.pdf'

%save('./ScanResults_uniform2V_wider');
