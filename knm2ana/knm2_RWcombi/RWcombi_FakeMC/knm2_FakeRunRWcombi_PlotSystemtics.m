
Mode = '3P';%''Gauss';

switch Mode
    case '3P'
        SysErr = [0.04:0.02:0.16,0.2];
        legstr = 'Step-wise plasma potential time evolution';
    case 'Gauss'
        SysErr = [0.01,0.02:0.02:0.1];
        legstr = 'Gaussian plasma potential time evolution';
end

%% load data
savedir  = [getenv('SamakPath'),'knm2ana/knm2_RWcombi/results/'];

savename = arrayfun(@(x) ...
    [savedir,sprintf('knm2_FakeRunRWcombi_Systematics_%s_%.0fmeVErr.mat',Mode,x)],SysErr*1e3,...
    'UniformOutput',0);
d = cellfun(@(x) importdata(x),savename,'UniformOutput',0);
%%
LmNuSq_sys = cell2mat(cellfun(@(x) sqrt(x.LmNuSq_cm^2-x.LmNuSq_stat^2),d,'UniformOutput',0));
LmNuSq_cm = cell2mat(cellfun(@(x) x.LmNuSq_cm,d,'UniformOutput',0));

%% plot
figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
plot(SysErr*1e3,LmNuSq_cm*1e3,'-o','LineWidth',2,'MarkerFaceColor',rgb('DodgerBlue'));
xlabel('Systematic plasma time evolution uncertainty (meV)');
ylabel(sprintf('1\\sigma sensitivity on (m_\\nu^2) (meV^2)'));
PrettyFigureFormat('FontSize',24);
xlim([min(SysErr)-0.002,max(SysErr)+0.002]*1e3)
leg = legend(legstr,'Location','northwest');
leg.EdgeColor = rgb('Silver');