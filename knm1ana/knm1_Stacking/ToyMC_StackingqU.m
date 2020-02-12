% choose distribution of T2 and DT
dist = 'gaussian';%'uniform';%'gaussian';
nRuns = 10000;

%% simulate single runs
switch dist
    case 'gaussian'
        qU   = (randn(nRuns,1).*0.001+1).*(18575-40);  % t2 concentrations: gaussian with mean=10 and std=2
    case 'uniform'
        qU   = rand(nRuns,1).*2+10;
    case '2gaussians'
        qU   = [randn(nRuns/2,1).*2+10;randn(nRuns/2,1).*3+30];
end
%Time = (1+randn(nRuns,1).*0.5).*(1919648/nRuns);    % time: gaussian 
Time = [(1+randn(nRuns/2,1).*0.005).*0.7.*(1919648/nRuns);(1+randn(nRuns/2,1).*0.005).*1.3.*(1919648/nRuns)];    % time: gaussian 
SingleRunCounts = ComputeCounts(qU,Time);
% %%
% fig1 = figure('Renderer','opengl');
% set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 0.5]);
%           
% h1 = histogram(qU);
% h1.FaceColor =rgb('DarkGreen');
% xlabel('qU (eV)');
% PrettyFigureFormat;
% save_name = [getenv('SamakPath'),'knm1ana/knm1Twins/plots/MC_qUDistribution_',dist,'.png'];
%print(fig1,save_name,'-dpng','-r400');
%% compute stacked run:
TimeStack = sum(Time);
%qUStack   = wmean(qU,Time);
%qUStack   = (sum(qU.^(-0.5).*Time)/sum(Time))^(0.5);

CountsStackData  = sum(SingleRunCounts);

qUStack =interp1(SingleRunCounts./Time,qU,CountsStackData./TimeStack,'spline');

CountsStackModel = ComputeCounts(qUStack,TimeStack);
%% compare
fprintf(2,'norm. residuals (stackdata - stackmodel) = %.10g sigma \n',(CountsStackData-CountsStackModel)./sqrt(CountsStackData))
%%
fig1 = figure('Renderer','opengl');
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.5, 0.5]);
psingle = plot(SingleRunCounts./Time,qU,'d','Color',rgb('DarkOliveGreen'),'MarkerFaceColor',rgb('DarkOliveGreen'));
hold on;
 plot(linspace(min(SingleRunCounts./Time),CountsStackData./TimeStack,numel(Time)),...
    qUStack.*ones(numel(Time),1),'--','LineWidth',3,'Color',rgb('GoldenRod'))
plot(CountsStackData./TimeStack.*ones(numel(Time),1),linspace(min(qU),qUStack,numel(Time)),...
    '--','LineWidth',2,'Color',rgb('GoldenRod'));
pstack = plot(CountsStackData./TimeStack,qUStack,'o','Color',rgb('Orange'),'MarkerFaceColor',rgb('GoldenRod'),...
    'MarkerSize',10);
%pmean = plot(CountsStackData./TimeStack,wmean(qU,Time),'s','Color',rgb('RoyalBlue'),'MarkerFaceColor',rgb('CadetBlue'),...
 %   'MarkerSize',10);
hold off;
PrettyFigureFormat;
xlabel('count rate (cps)');
ylabel('qU (eV)');
ylim([min(qU),max(qU)]);
xlim([min(SingleRunCounts./Time),max(SingleRunCounts./Time)]);
leg = legend([psingle,pstack],'single runs','stacked runs');%,'weighted average');
legend boxoff
set(gca,'FontSize',18);
 save_name = [getenv('SamakPath'),'knm1ana/knm1Twins/plots/MC_qUInverted','.png'];
%print(fig1,save_name,'-dpng','-r400');
%close
%%
function Counts = ComputeCounts(qu,time) 
%Counts = arrayfun(@(x,t) x.^(-1).*t.*2.3*1e3,qu,time);
Counts = arrayfun(@(x,t) (-x.*0.5+4e5).*t.*2.3*1e3,qu,time);
end
