% Get std of Bana and qU for kNM2 golden runs

savedir = [getenv('SamakPath'),'knm2ana/knm2_UniformFit/results/'];
savename = sprintf('%sknm2_AnalyzingPlaneInhomo.mat',savedir);

if exist(savename,'file')
    load(savename)
else
RunAnaArg = {'RunList','KNM2_Prompt',... % all KNM2 golden runs
    'fixPar','mNu E0 Bkg Norm',...           % free Parameter !!
    'DataType','Real',...
    'AnaFlag','StackPixel'};

M = MultiRunAnalysis(RunAnaArg{:});
M.ReadSingleRunData;
RunList = M.RunList;
%% retarding potential:
% this value from from field-maps -> the same for all runs and subruns
qU    = M.SingleRunData.qU;
qUStd =  squeeze(mean(mean(std(qU,0,3)))); % standard deviation over pixels

%% B-ana: same for all runs
Ba = M.SingleRunData.MACE_Ba_T;
BaStd = mean(std(Ba));

%% display
fprintf('KNM2 analyzing plane inhomogeneities----------------------------\n')
fprintf('qU: sigma = %.0f meV    --> m-nu approx. %.0e eV^2 \n',qUStd*1e3,-2*qUStd^2);
fprintf('Ba: sigma = %.1e T  \n',BaStd);
fprintf('----------------------------------------------------------------\n')


%% plot broadened response function
RFdef = M.ModelObj.RF;
Te    = M.ModelObj.Te;
qUav    = M.ModelObj.qU;

M.ModelObj.MACE_Sigma = qUStd;
M.ModelObj.InitializeRF;
RFbro = M.ModelObj.RF;

MakeDir(savedir);
save(savename,'Te','qU','RFdef','RFbro','qUStd','qUav','Ba','BaStd','RunAnaArg','RunList');
end
%% plot
f2 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
qUi = 20;
pref = plot(Te-qUav(qUi),RFdef(:,qUi),'-','LineWidth',2,'Color',rgb('DodgerBlue'));
hold on;
pbro =  plot(Te-qUav(qUi),RFbro(:,qUi),'-.','LineWidth',2,'Color',rgb('Orange'));
leg = legend([pref,pbro],sprintf('\\sigma = 0 meV'),sprintf('\\sigma = %.0f meV',qUStd*1e3),...
    'EdgeColor',rgb('Silver'),'Location','northwest');
xlabel('Surplus energy (eV)');
ylabel('Transmission probability');
PrettyFigureFormat('FontSize',24);
xlim([-1 40]);
plotname = strrep(strrep(savename,'results','plots'),'.mat','.pdf');
export_fig(f2,plotname);

xlim([2.7 3.2]);
export_fig(f2,strrep(plotname,'.pdf','_Zoom.pdf'));

