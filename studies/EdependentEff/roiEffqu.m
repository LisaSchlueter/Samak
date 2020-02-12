%
% Energy dependent efficiency
% Data by Sanshiro / L�o
%
% ROIqU: Retarding Potential
% ROIEff: Efficiency Correction
%
% 
% Last Updated: 26/07/2018
%

format long

% Retarding Potential
ROIqU = [14.6        14.61        14.62        14.63        14.64        14.65        14.66        14.67        14.68        14.69         14.7        14.71        14.72        14.73        14.74        14.75        14.76        14.77        14.78        14.79         14.8        14.81        14.82        14.83        14.84        14.85        14.86        14.87        14.88        14.89         14.9        14.91        14.92        14.93        14.94        14.95        14.96        14.97        14.98        14.99           15        15.01        15.02        15.03        15.04        15.05        15.06        15.07        15.08        15.09         15.1        15.11        15.12        15.13        15.14        15.15        15.16        15.17        15.18        15.19         15.2        15.21        15.22        15.23        15.24        15.25        15.26        15.27        15.28        15.29         15.3        15.31        15.32        15.33        15.34        15.35        15.36        15.37        15.38        15.39         15.4        15.41        15.42        15.43        15.44        15.45        15.46        15.47        15.48        15.49         15.5        15.51        15.52        15.53        15.54        15.55        15.56        15.57        15.58        15.59         15.6        15.61        15.62        15.63        15.64        15.65        15.66        15.67        15.68        15.69         15.7        15.71        15.72        15.73        15.74        15.75        15.76        15.77        15.78        15.79         15.8        15.81        15.82        15.83        15.84        15.85        15.86        15.87        15.88        15.89         15.9        15.91        15.92        15.93        15.94        15.95        15.96        15.97        15.98        15.99           16        16.01        16.02        16.03        16.04        16.05        16.06        16.07        16.08        16.09         16.1        16.11        16.12        16.13        16.14        16.15        16.16        16.17        16.18        16.19         16.2        16.21        16.22        16.23        16.24        16.25        16.26        16.27        16.28        16.29         16.3        16.31        16.32        16.33        16.34        16.35        16.36        16.37        16.38        16.39         16.4        16.41        16.42        16.43        16.44        16.45        16.46        16.47        16.48        16.49         16.5        16.51        16.52        16.53        16.54        16.55        16.56        16.57        16.58        16.59         16.6        16.61        16.62        16.63        16.64        16.65        16.66        16.67        16.68        16.69         16.7        16.71        16.72        16.73        16.74        16.75        16.76        16.77        16.78        16.79         16.8        16.81        16.82        16.83        16.84        16.85        16.86        16.87        16.88        16.89         16.9        16.91        16.92        16.93        16.94        16.95        16.96        16.97        16.98        16.99           17        17.01        17.02        17.03        17.04        17.05        17.06        17.07        17.08        17.09         17.1        17.11        17.12        17.13        17.14        17.15        17.16        17.17        17.18        17.19         17.2        17.21        17.22        17.23        17.24        17.25        17.26        17.27        17.28        17.29         17.3        17.31        17.32        17.33        17.34        17.35        17.36        17.37        17.38        17.39         17.4        17.41        17.42        17.43        17.44        17.45        17.46        17.47        17.48        17.49         17.5        17.51        17.52        17.53        17.54        17.55        17.56        17.57        17.58        17.59         17.6        17.61        17.62        17.63        17.64        17.65        17.66        17.67        17.68        17.69         17.7        17.71        17.72        17.73        17.74        17.75        17.76        17.77        17.78        17.79         17.8        17.81        17.82        17.83        17.84        17.85        17.86        17.87        17.88        17.89         17.9        17.91        17.92        17.93        17.94        17.95        17.96        17.97        17.98        17.99           18        18.01        18.02        18.03        18.04        18.05        18.06        18.07        18.08        18.09         18.1        18.11        18.12        18.13        18.14        18.15        18.16        18.17        18.18        18.19         18.2        18.21        18.22        18.23        18.24        18.25        18.26        18.27        18.28        18.29         18.3        18.31        18.32        18.33        18.34        18.35        18.36        18.37        18.38        18.39         18.4        18.41        18.42        18.43        18.44        18.45        18.46        18.47        18.48        18.49         18.5        18.51        18.52        18.53        18.54        18.55        18.56        18.57        18.58        18.59         18.6        18.61];
[i,j]=find(ROIqU>16);
ROIqU = ROIqU(j);

% Energy dependent efficiency Correction
ROIEff = [0.98833     0.98833      0.9884     0.98844     0.98848     0.98851     0.98851     0.98858     0.98858     0.98866     0.98869     0.98873     0.98876      0.9888     0.98884     0.98884     0.98891     0.98894     0.98898     0.98901     0.98905     0.98909     0.98908     0.98916     0.98919     0.98923     0.98923      0.9893     0.98933     0.98937      0.9894      0.9894     0.98947     0.98947     0.98955     0.98958     0.98961     0.98965     0.98968     0.98972     0.98975     0.98979     0.98982     0.98986     0.98989     0.98993     0.98996     0.98999     0.99003     0.99006      0.9901     0.99013     0.99017      0.9902     0.99024     0.99027      0.9903     0.99034     0.99037      0.9904     0.99044     0.99047     0.99051     0.99054     0.99057     0.99061     0.99064     0.99068     0.99071     0.99074     0.99077     0.99081     0.99084     0.99088     0.99091     0.99094     0.99098     0.99101     0.99104     0.99108     0.99111     0.99114     0.99118     0.99121     0.99124     0.99127     0.99131     0.99134     0.99137     0.99141     0.99144     0.99147      0.9915     0.99154     0.99157      0.9916     0.99164     0.99167      0.9917     0.99173     0.99177      0.9918     0.99183     0.99186      0.9919     0.99193     0.99196     0.99199     0.99203     0.99206     0.99209     0.99212     0.99216     0.99219     0.99222     0.99225     0.99228     0.99232     0.99235     0.99238     0.99241     0.99244     0.99248     0.99251     0.99254     0.99257      0.9926     0.99264     0.99267      0.9927     0.99273     0.99276     0.99279     0.99283     0.99286     0.99289     0.99292     0.99295     0.99298     0.99301     0.99305     0.99308     0.99311     0.99314     0.99317     0.99321     0.99324     0.99327      0.9933     0.99333     0.99336     0.99339     0.99343     0.99346     0.99349     0.99352     0.99355     0.99358     0.99361     0.99365     0.99368     0.99371     0.99374     0.99377      0.9938     0.99383     0.99386     0.99389     0.99392     0.99396     0.99399     0.99402     0.99405     0.99408     0.99411     0.99414     0.99417      0.9942     0.99423     0.99426      0.9943     0.99433     0.99436     0.99439     0.99442     0.99445     0.99448     0.99451     0.99454     0.99457      0.9946     0.99463     0.99466     0.99469     0.99472     0.99476     0.99479     0.99482     0.99485     0.99488     0.99491     0.99494     0.99497       0.995     0.99503     0.99506     0.99509     0.99512     0.99515     0.99518     0.99521     0.99524     0.99527      0.9953     0.99533     0.99536     0.99539     0.99542     0.99545     0.99548     0.99551     0.99554     0.99557      0.9956     0.99563     0.99567     0.99569     0.99572     0.99575     0.99578     0.99581     0.99584     0.99587      0.9959     0.99593     0.99596     0.99599     0.99602     0.99605     0.99608     0.99611     0.99614     0.99617      0.9962     0.99623     0.99626     0.99629     0.99632     0.99635     0.99637      0.9964     0.99643     0.99646     0.99649     0.99652     0.99655     0.99658     0.99661     0.99664     0.99667      0.9967     0.99673     0.99676     0.99679     0.99682     0.99685     0.99688     0.99691     0.99694     0.99697       0.997     0.99703     0.99706     0.99708     0.99711     0.99714     0.99717      0.9972     0.99723     0.99726     0.99729     0.99732     0.99735     0.99738     0.99741     0.99744     0.99746     0.99749     0.99752     0.99755     0.99758     0.99761     0.99764     0.99767      0.9977     0.99772     0.99775     0.99778     0.99781     0.99784     0.99787      0.9979     0.99792     0.99795     0.99798     0.99801     0.99804     0.99806     0.99809     0.99812     0.99815     0.99818     0.99821     0.99823     0.99826     0.99828     0.99831     0.99834     0.99837      0.9984     0.99842     0.99845     0.99848      0.9985     0.99853     0.99855     0.99858     0.99861     0.99863     0.99866     0.99869     0.99871     0.99874     0.99877     0.99879     0.99882     0.99884     0.99887     0.99889     0.99892     0.99894     0.99897     0.99899     0.99902     0.99905     0.99907     0.99909     0.99912     0.99914     0.99916     0.99919     0.99921     0.99923     0.99925     0.99928      0.9993     0.99932     0.99935     0.99937     0.99939     0.99941     0.99943     0.99945     0.99947     0.99949     0.99952     0.99953     0.99955     0.99957      0.9996     0.99961     0.99963     0.99965     0.99967     0.99968      0.9997     0.99973     0.99974     0.99975     0.99977     0.99978      0.9998     0.99981     0.99984     0.99984     0.99985     0.99986     0.99987     0.99988     0.99989     0.99992     0.99991     0.99992     0.99993     0.99996     0.99994     0.99995     0.99996     0.99996     0.99997     0.99997           1]';
ROIEff = ROIEff(j);


% Interpolation
ROIEffInterp = @(qu) interp1(ROIqU,ROIEff,qu,'spline');
qu = linspace(min(ROIqU),max(ROIqU),1000);

% Plot
p=Plot(ROIqU,ROIEff,qu,ROIEffInterp(qu));
p.XLabel = 'qU (V)';
p.YLabel = 'Correction Factor';
p.Title  = 'Energy dependent Efficiency';

%% Fit: 'qU  Dependent FPD Efficiency'.
[xData, yData] = prepareCurveData( ROIqU, ROIEff );

% Set up fittype and options.
ft = fittype( 'poly4' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure(1)
errorbar(ROIqU,ROIEff,ROIEff*0,'-s','MarkerSize',5,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue');
hold on
plot(ROIqU,fitresult(ROIqU),'LineWidth',2);
hold off
xlabel('qU (kV)');
ylabel('Correction Factor');
legend('data','4-order polynomial fit','Location','NorthWest')
title('ROI - qU Dependent Efficiency - Fit');
grid on
PrettyFigureFormat
