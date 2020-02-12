% choose distribution of T2 and DT
dist = '2gaussians';%'uniform';%'gaussian';
nRuns = 10000;

%% simulate single runs
switch dist
    case 'gaussian'
        T2   = (1+randn(nRuns,1).*0.05).*0.95;  % t2 concentrations: gaussian with mean=10 and std=2
        DT   = (1+randn(nRuns,1).*0.05).*0.01; 
    case 'uniform'
        T2   = rand(nRuns,1).*2*0.05-0.05+0.95;
        DT   = rand(nRuns,1).*2*5e-04-5e-04+0.01;
    case '2gaussians'
      T2   = [(1+randn(nRuns/2,1).*0.01).*0.95;(1+randn(nRuns/2,1).*0.01).*0.91];
      DT   = [(1+randn(nRuns/2,1).*0.01).*0.0098;(1+randn(nRuns/2,1).*0.01).*0.0102];
end
Time = randn(nRuns,1)+1919648/nRuns;    % time: gaussian with mean 100
SingleRunCounts = ComputeCounts(T2,DT,Time);
%%
fig1 = figure('Renderer','opengl');
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 1, 0.5]);
           subplot(1,2,1);
h1 = histogram(T2);
h1.FaceColor =rgb('DarkGreen');
xlabel('mol. frac. T_2');
PrettyFigureFormat;

subplot(1,2,2);
h1 = histogram(DT);
h1.FaceColor =rgb('DarkGreen');
xlabel('mol. frac. DT');
PrettyFigureFormat;

save_name = [getenv('SamakPath'),'knm1ana/knm1Twins/plots/MC_IsotopologueDistribution_',dist,'.png'];
%print(fig1,save_name,'-dpng','-r400');
%% compute stacked run:
TimeStack = sum(Time);
T2Stack   = wmean(T2,Time);
DTStack   = wmean(DT,Time);
CountsStackData  = sum(SingleRunCounts);
CountsStackModel = ComputeCounts(T2Stack,DTStack,TimeStack);

%% compare
fprintf(2,'norm. residuals (stackdata - stackmodel) = %.10g sigma \n',(CountsStackData-CountsStackModel)./sqrt(CountsStackData))

function Counts = ComputeCounts(t2,dt,time) 
Counts = arrayfun(@(x,y,z) (x+0.5.*y).*z,t2,dt,time);
end
%% since all stacked SC increase the number of counts linearly over TIME
% stacking with respect to time is not an issue!