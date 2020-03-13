% plot final states and error band

% config
DataType ='Real';
RunList = 'KNM1';
exclDataStart = 13;
chi2 = 'chi2Stat';

% model
CommongArg = {'RunList',RunList,'chi2','chi2Stat','DataType',DataType,...
    'fixPar','5 6 7 8 9 10 11','exclDataStart',exclDataStart,'chi2',chi2,...
    'NonPoissonScaleFactor',1};
if ~exist('M','var')
    M = MultiRunAnalysis(CommongArg{:});
end

M.ComputeCM('SysEffect',struct('FSD','ON'),'BkgCM','OFF','nTrials',5000,'RecomputeFlag','ON');

CM = M.FitCM_Obj; CM.RecomputeFlag = 'OFF';
CM.ComputeCM_FSD;

d = importdata(CM.CovMatFile);
%%
Isotopologue = 'TT';
FSD_P = d.([Isotopologue,'_P_norm'])'; %probabilities
FSD_E = d.obj.StudyObject.([Isotopologue,'exE']);

f22 = figure('Renderer','opengl');
set(f22, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.5]);
[l, a] = boundedline(FSD_E,mean(FSD_P),std(FSD_P));
l.Color = rgb('FireBrick'); l.LineWidth = 2;
a.FaceColor =rgb('FireBrick'); a.FaceAlpha = 0.3;
PrettyFigureFormat;
xlabel('excitation energy (eV)');
ylabel('probability');
xlim([0,90]);
set(gca,'XScale','log');
switch Isotopologue
    case 'TT'
        is_str = sprintf('T_2');
end
leg = legend([l,a],sprintf('%s - %s',is_str,M.ModelObj.([Isotopologue,'FSD'])),sprintf('1\\sigma error band'));
legend boxoff;

savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/plots/'];
MakeDir(savedir);
print(f22,[savedir,'FSDsystematics.png'],'-dpng','-r450');


