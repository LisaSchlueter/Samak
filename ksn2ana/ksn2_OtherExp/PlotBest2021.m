% plot result from BEST experiment
% arxiv: arxiv.org/abs/2109.11482 
% extracted with DataThief software

% load data from file
dirBest = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm2/Others/'];
fBest1a = [dirBest,'contour_Best_1sigma_a.txt'];
fBest1b = [dirBest,'contour_Best_1sigma_b.txt'];
fBest2a = [dirBest,'contour_Best_2sigma_a.txt'];
fBest2b = [dirBest,'contour_Best_2sigma_b.txt'];
fBest3  = [dirBest,'contour_Best_3sigma.txt'];

dataBest1a = importdata(fBest1a);
dataBest1b = importdata(fBest1b);
dataBest2a = importdata(fBest2a);
dataBest2b = importdata(fBest2b);
dataBest3 = importdata(fBest3);

% plot
GetFigure;
contour_1a = plot(dataBest1a(:,1),dataBest1a(:,2),'-r');
hold on;
contour_1b = plot(dataBest1b(:,1),dataBest1b(:,2),'-r');
contour_2a = plot(dataBest2a(:,1),dataBest2a(:,2),'-b');
contour_2b = plot(dataBest2b(:,1),dataBest2b(:,2),'-b');
contour_3 = plot(dataBest3(:,1),dataBest3(:,2),'-g');

xlabel(sprintf('sin^2(2\\theta_{ee})'));
ylabel(sprintf('\\Deltam^2 (eV^2)'));
PrettyFigureFormat;
leg = legend([contour_1a,contour_2a,contour_3],sprintf('1\\sigma'),sprintf('2\\sigma'),sprintf('3\\sigma'));
PrettyLegendFormat(leg);
