% fit  with different final state distributions
% plot impact on neutrino mass

% config
DataType ='Twin';
RunList = 'KNM1';
exclDataStart = 14;
chi2 = 'chi2Stat';
ReFit = 'OFF';

%label results and load if possible
savedir = [getenv('SamakPath'),'knm1ana/knm1_systematics/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('NuMassShift_FSD_%s_%s_%s_%.0f.mat',RunList,DataType,chi2,exclDataStart)];

if exist(savename,'file') && strcmp(ReFit,'OFF')
    load(savename)
else
CommongArg = {'RunList',RunList,'chi2','chi2Stat','DataType',DataType,...
    'fixPar','5 6 7 8 9 10 11','exclDataStart',exclDataStart,'chi2',chi2,...
    'minuitOpt','min;minos'};

M = MultiRunAnalysis(CommongArg{:});

mNuSq = zeros(5,1); mNuSqErr = zeros(5,1);

M.ModelObj.TTFSD = 'Sibille';
M.ModelObj.LoadFSD;
M.Fit;
mNuSq(1)    = M.FitResult.par(1); 
mNuSqErr(1) = M.FitResult.err(1);

M.ModelObj.TTFSD = 'Sibille0p5eV';
M.ModelObj.LoadFSD;
M.Fit;
mNuSq(2)    = M.FitResult.par(1); 
mNuSqErr(2) = M.FitResult.err(1);

M.ModelObj.TTFSD = 'SAENZ';
M.ModelObj.LoadFSD;
M.Fit;
mNuSq(3)    = M.FitResult.par(1); 
mNuSqErr(3) = M.FitResult.err(1);

M.ModelObj.TTFSD = 'HT';
M.ModelObj.LoadFSD;
M.Fit;
mNuSq(4)    = M.FitResult.par(1); 
mNuSqErr(4) = M.FitResult.err(1);

M.ModelObj.TTFSD = 'DOSS';
M.ModelObj.LoadFSD;
M.Fit;
mNuSq(5)    = M.FitResult.par(1); 
mNuSqErr(5) = M.FitResult.err(1);

M.ModelObj.TTFSD = 'ROLL';
M.ModelObj.LoadFSD;
M.Fit;
mNuSq(6)     = M.FitResult.par(1); 
mNuSqErr(6) = M.FitResult.err(1);

FSDname = {'Sibille','Sibille0p5','Saenz', 'Saenz HT','Doss','Roll'};
save(savename,'mNuSq','mNuSqErr','FSDname');
end

%% plot neutrino mass shift as a function of FSD
f44 = figure('Renderer','opengl');
set(f44, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.5]);

select = 1:6;
y    = mNuSq(select);
yErr = mNuSqErr(select);
xstr = {'Sibille','Sibille rebinned','Saenz', 'Doss','Roll','Saenz HT'};
ystr = sprintf('m_\\nu^2 - < m_\\nu^2 > (eV^2)');

x = 1:numel(y);

p =  plot(x,y-mean(y),'--o','LineWidth',3,'MarkerSize',8,'Color',M.PlotColor);
ylabel(ystr);
xticks(x);
xticklabels(xstr);
PrettyFigureFormat;
xlabel('final state distribution')
set(gca,'FontSize',24);
grid on
plotname = strrep(savename,'results','plots');
print(f44,plotname,'-dpng','-r450');
