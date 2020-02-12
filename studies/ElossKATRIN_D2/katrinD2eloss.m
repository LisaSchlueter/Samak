% Compute and Display KATRIN ELoss Function
% Th. Lasserre, 12/03/2019

% ELoss Function Parameters
amp1 	= 0.0288017; 
pos1 	= 11.8468; 
sig1 	= 0.18288;
amp2 	= 0.266339;
pos2 	= 12.785; 	
sig2 	= 0.486158; 
amp3 	= 0.0728571;
pos3 	= 14.9051; 	
sig3 	= 1.0307; 	
amptof 	= 1; 		
mutof 	= 0.104427; 
norm 	= 0.997937; 
mu1 	= 0.0964532;
mu2 	= 0.656861; 
mu3 	= 1.64586; 	
tail 	= 1; 	
n1 	    = 1.0099;

% Covariance Matrix of Fit Parameters
elossKATRINcm = ...
[ 3.13e-06  1.06e-05  3.80e-06 -2.28e-08  2.76e-06 -4.21e-06 -8.42e-08 -3.75e-06  1.10e-05  3.29e-10 -1.48e-09  4.88e-10  1.92e-09  4.09e-09  7.47e-11;
 1.06e-05  3.07e-04  2.01e-04 -4.59e-06  3.77e-05 -2.72e-05 -3.81e-07 -1.31e-05  4.33e-05  2.37e-11  4.05e-08 -1.18e-09 -1.31e-08 -8.22e-08  2.39e-10;
 3.80e-06  2.01e-04  1.92e-04 -1.84e-06  2.90e-05 -2.76e-05 -2.31e-07 -1.84e-05  4.51e-05  1.29e-09  1.55e-08  1.22e-09  1.90e-09 -1.69e-08  3.60e-10;
-2.28e-08 -4.59e-06 -1.84e-06  1.57e-06  1.03e-06  5.77e-09  3.19e-07  2.76e-06 -2.07e-05 -1.78e-08 -2.91e-07  1.08e-08 -1.04e-08 -1.80e-09  9.32e-09;
 2.76e-06  3.77e-05  2.90e-05  1.03e-06  1.48e-05 -1.43e-06  3.75e-07  1.08e-05 -3.81e-05 -1.81e-10  2.74e-08 -5.93e-09 -1.48e-08 -5.69e-09 -2.68e-09;
-4.21e-06 -2.72e-05 -2.76e-05  5.77e-09 -1.43e-06  1.56e-05  5.77e-07  1.98e-05 -6.46e-05 -7.16e-09 -1.43e-07  5.57e-09 -1.86e-09  1.43e-08  4.52e-09;
-8.42e-08 -3.81e-07 -2.31e-07  3.19e-07  3.75e-07  5.77e-07  2.26e-07 -1.92e-07 -1.08e-05 -6.42e-09 -1.21e-07  4.30e-09 -4.88e-09  1.52e-10  3.82e-09;
-3.75e-06 -1.31e-05 -1.84e-05  2.76e-06  1.08e-05  1.98e-05 -1.92e-07  6.59e-05 -1.07e-04  5.48e-08  1.12e-06 -4.31e-08  3.79e-08 -2.81e-08 -3.70e-08;
 1.10e-05  4.33e-05  4.51e-05 -2.07e-05 -3.81e-05 -6.46e-05 -1.08e-05 -1.07e-04  7.30e-04  1.70e-07  3.59e-06 -1.38e-07  1.17e-07 -8.35e-08 -1.18e-07;
 3.29e-10  2.37e-11  1.29e-09 -1.78e-08 -1.81e-10 -7.16e-09 -6.42e-09  5.48e-08  1.70e-07  6.44e-08  5.06e-08 -1.58e-09  2.40e-09  1.90e-09 -1.46e-09;
-1.48e-09  4.05e-08  1.55e-08 -2.91e-07  2.74e-08 -1.43e-07 -1.21e-07  1.12e-06  3.59e-06  5.06e-08  7.37e-07 -2.51e-08  3.44e-08  2.72e-08 -2.27e-08;
 4.88e-10 -1.18e-09  1.22e-09  1.08e-08 -5.93e-09  5.57e-09  4.30e-09 -4.31e-08 -1.38e-07 -1.58e-09 -2.51e-08  2.09e-08 -1.16e-09 -1.04e-09  9.05e-09;
 1.92e-09 -1.31e-08  1.90e-09 -1.04e-08 -1.48e-08 -1.86e-09 -4.88e-09  3.79e-08  1.17e-07  2.40e-09  3.44e-08 -1.16e-09  3.83e-08  4.76e-09 -1.18e-09;
 4.09e-09 -8.22e-08 -1.69e-08 -1.80e-09 -5.69e-09  1.43e-08  1.52e-10 -2.81e-08 -8.35e-08  1.90e-09  2.72e-08 -1.04e-09  4.76e-09  1.86e-08 -1.25e-09;
 7.47e-11  2.39e-10  3.60e-10  9.32e-09 -2.68e-09  4.52e-09  3.82e-09 -3.70e-08 -1.18e-07 -1.46e-09 -2.27e-08  9.05e-09 -1.18e-09 -1.25e-09  6.62e-09;];

% Compute ELoss
e        = [0.:.1:60]; % eV
f        = en_loss_sc(e, amp1, pos1, sig1, amp2, pos2, sig2, amp3, pos3, sig3);

% Generate ELoss Envelop
pv        = [amp1, pos1, sig1, amp2, pos2, sig2, amp3, pos3, sig3, mutof, norm, mu1, mu2, mu3, n1];
mpv       = mvnrnd(pv,elossKATRINcm,10000);
mf        = zeros(10000,numel(e));
for i=1:1:10000
mf(i,:) = en_loss_sc(e, mpv(i,1), mpv(i,2), mpv(i,3), mpv(i,4), mpv(i,5), mpv(i,6), mpv(i,7), mpv(i,8), mpv(i,9));
end
fe      = cov(mf);

%% Plot Function & Error
fig1 = figure('Renderer','openGL');
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.9, 0.9]);
s1=subplot(2,2,1)
[lb1 , pb1]  = boundedline(e,f,sqrt(diag(fe))*10,'alpha','transparency', 0.4,'cmap',rgb('Amethyst')); 
lb1.LineStyle= '--';
hold on;
plot(e,f,'LineWidth',2,'Color',rgb('SteelBlue'));
hold off
xlabel('E-qU (eV)');
ylabel('energy loss function');
PrettyFigureFormat
s2=subplot(2,2,3)
[lb2 , pb2]  = boundedline(e,zeros(1,numel(e)),(sqrt(diag(fe))'),'alpha','transparency', 0.4,'cmap',rgb('Amethyst')); 
xlabel('E-qU (eV)');
ylabel('confidence interval');
PrettyFigureFormat
subplot(2,2,[2 4])
corplot(elossKATRINcm)
title('KATRIN energy loss function - D_2')
linkaxes([s1,s2],'x');

