% systematic breakdown for magnetic fields
exclDataStart = 14;
DataType = 'Twin';
savedir = [getenv('SamakPath'),'knm1ana/knm1sensitivity/results/'];
savename = [savedir,sprintf('knm1_sensitivityBFBreakdown_%s_%.0f.mat',DataType,exclDataStart)];

if exist(savename,'file')
    load(savename)
else
M = MultiRunAnalysis('RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'FSDFlag','Sibille0p5eV','exclDataStart',exclDataStart,'chi2','chi2Stat','fitter','matlab');
M.chi2 = 'chi2Stat';
M.Fit;
mNuSqErrStat = M.FitResult.err(1);
mNuSqErrCM = zeros(4,1);

M.chi2 = 'chi2CMShape';
% all bf
M.ComputeCM('SysEffect',struct('RF_BF','ON'),...
'nTrials',1000,'BkgCM','OFF');
M.Fit;
mNuSqErrCM(4) = M.FitResult.err(1);

%bmax (pinch)
M.ComputeCM('SysEffect',struct('RF_BF','ON'),...
'MACE_Bmax_T_RelErr',0.002,...
'WGTS_B_T_RelErr',0,...
'MACE_Ba_T_RelErr',0,...
'nTrials',1000,'BkgCM','OFF');
M.Fit;
mNuSqErrCM(1) = M.FitResult.err(1);

% b source
M.ComputeCM('SysEffect',struct('RF_BF','ON'),...
'MACE_Bmax_T_RelErr',0,...
'WGTS_B_T_RelErr',0.025,... 
'MACE_Ba_T_RelErr',0,...
'nTrials',1000,'BkgCM','OFF');
M.Fit;
mNuSqErrCM(2) = M.FitResult.err(1);


% ba (analyzing plane)
M.ComputeCM('SysEffect',struct('RF_BF','ON'),...
'MACE_Bmax_T_RelErr',0,...
'WGTS_B_T_RelErr',0,... 
'MACE_Ba_T_RelErr',0.01,...
'nTrials',1000,'BkgCM','OFF');
M.Fit;
mNuSqErrCM(3) = M.FitResult.err(1);

save(savename,'mNuSqErrCM','mNuSqErrStat');
end
%% get sys contribution & plot
mNuSqErrSys = real(sqrt(mNuSqErrCM.^2-mNuSqErrStat.^2));
SingleBarX = linspace(5.393,10.85,numel(mNuSqErrSys));
bsingle = cell(numel(mNuSqErrCM),1);

f33 = figure('Renderer','opengl');
set(f33, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);

PlotColor = {'DodgerBlue','FireBrick','GoldenRod','Silver'};
for i=1:numel(mNuSqErrSys)
    bsingle{i}  = barh(SingleBarX(i)',mNuSqErrSys(i)','stacked','BarWidth',1.5);
    hold on;
    btmp= bsingle{i};
    btmp.FaceColor = rgb(PlotColor{i}); btmp.LineStyle ='none';
    btmp.FaceAlpha=0.8;
end

hold off;
PrettyFigureFormat;
yticks(SingleBarX);
grid on;
yticklabels({sprintf('B_{max}'),sprintf('B_s'),sprintf('B_a'),'combined'})
xlabel(sprintf('m_{\\beta}^2 sensitivity at 68.3 %% C.L.'));
title('magnetic fields - systematics breakdown');

plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
publish_figurePDF(f33,plotname);

