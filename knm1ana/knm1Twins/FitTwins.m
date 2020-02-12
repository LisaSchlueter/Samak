function [mNuSq,mNuSqErr,E0,E0Err,N,NErr,B,BErr,dof,chi2min] = FitTwins(varargin)

% function to fit all twins (with different same flags) and save fit results

p=inputParser;
p.addParameter('RunList','KNM1',@(x)ischar(x));
p.addParameter('exclDataStart',14,@(x)isfloat(x));
p.addParameter('ReFit','OFF',@(x)ischar(x));
p.addParameter('SavePlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.parse(varargin{:});

RunList       = p.Results.RunList;
exclDataStart = p.Results.exclDataStart;
ReFit         = p.Results.ReFit;
SavePlot      = p.Results.SavePlot;

[Twin, TwinqU, TwinqUfrac, TwinqUqUfrac, TwinIs, TwinCD,Twinall] = ...
    ComputeLoadTwinObjects('RunList',RunList);

%% Fit
savedir = [getenv('SamakPath'),'knm1ana/knm1Twins/results/'];
FitRange = round(-Twin.ModelObj.qU(exclDataStart)+Twin.ModelObj.Q_i); %for labeling

save_file = [savedir,'FitResult_Twins_all',RunList,sprintf('%.0feVrange',FitRange),'.mat'];
if exist(save_file,'file') && strcmp(ReFit,'OFF')
    load(save_file,'mNuSq','mNuSqErr','E0','E0Err','N','NErr','B','BErr','dof','chi2min');
else
    
    % fit range
    Twin.exclDataStart = exclDataStart;
    TwinqU.exclDataStart = exclDataStart;
    TwinqUfrac.exclDataStart = exclDataStart;
    TwinCD.exclDataStart = exclDataStart;
    TwinIs.exclDataStart = exclDataStart;
    Twinall.exclDataStart = exclDataStart;
    TwinqUqUfrac.exclDataStart = exclDataStart;

    % fit parameter
    fixPar = '5 6 7 8 9 10 11'; % fix FSD parameter and qU-offset. free mNuSq, E0, B, N
    Twin.fixPar = fixPar;
    TwinqU.fixPar = fixPar;
    TwinqUfrac.fixPar = fixPar;
    TwinCD.fixPar = fixPar;
    TwinIs.fixPar = fixPar;
    Twinall.fixPar = fixPar;
    TwinqUqUfrac.fixPar = fixPar;

    % fit
    Twin.Fit;
    TwinqU.Fit;
    TwinqUfrac.Fit;
    TwinqUqUfrac.Fit;
    TwinCD.Fit;
    TwinIs.Fit;
    Twinall.Fit;

    
    mNuSq      = [Twin.FitResult.par(1),TwinCD.FitResult.par(1),TwinIs.FitResult.par(1),TwinqUfrac.FitResult.par(1),...
        TwinqU.FitResult.par(1),TwinqUqUfrac.FitResult.par(1),Twinall.FitResult.par(1)];%,TwinallbutqUfrac.FitResult.par(1)
    mNuSqErr   = [Twin.FitResult.err(1),TwinCD.FitResult.err(1),TwinIs.FitResult.err(1),TwinqUfrac.FitResult.err(1),...
        TwinqU.FitResult.err(1),TwinqUqUfrac.FitResult.err(1),Twinall.FitResult.err(1)];%,TwinallbutqUfrac.FitResult.err(1)
    
    E0         = [Twin.FitResult.par(2),TwinCD.FitResult.par(2),TwinIs.FitResult.par(2),TwinqUfrac.FitResult.par(2),...%TwinallbutqUfrac.FitResult.par(2),
        TwinqU.FitResult.par(2),TwinqUqUfrac.FitResult.par(2),Twinall.FitResult.par(2)];
    E0Err      = [Twin.FitResult.err(2),TwinCD.FitResult.err(2),TwinIs.FitResult.err(2),TwinqUfrac.FitResult.err(2),...
        TwinqU.FitResult.err(2),TwinqUqUfrac.FitResult.err(2),Twinall.FitResult.err(2)];%TwinallbutqUfrac.FitResult.err(2),
    
    B         = [Twin.FitResult.par(3),TwinCD.FitResult.par(3),TwinIs.FitResult.par(3),TwinqUfrac.FitResult.par(3),...
        TwinqU.FitResult.par(3),TwinqUqUfrac.FitResult.par(3),Twinall.FitResult.par(3)];%TwinallbutqUfrac.FitResult.par(3),
    BErr      = [Twin.FitResult.err(3),TwinCD.FitResult.err(3),TwinIs.FitResult.err(3),TwinqUfrac.FitResult.err(3),...
        TwinqU.FitResult.err(3),TwinqUqUfrac.FitResult.err(3),Twinall.FitResult.err(3)];%TwinallbutqUfrac.FitResult.err(3),
    
    N         = [Twin.FitResult.par(4),TwinCD.FitResult.par(4),TwinIs.FitResult.par(4),TwinqUfrac.FitResult.par(4),...
        TwinqU.FitResult.par(4),TwinqUqUfrac.FitResult.par(4),Twinall.FitResult.par(4)];%TwinallbutqUfrac.FitResult.par(4),
    NErr      = [Twin.FitResult.err(4),TwinCD.FitResult.err(4),TwinIs.FitResult.err(4),TwinqUfrac.FitResult.err(4),...
        TwinqU.FitResult.err(4),TwinqUqUfrac.FitResult.err(4),Twinall.FitResult.err(4)];%,TwinallbutqUfrac.FitResult.err(4)
    
    chi2min   = [Twin.FitResult.chi2min,TwinCD.FitResult.chi2min,TwinIs.FitResult.chi2min,TwinqUfrac.FitResult.chi2min,...
        TwinqU.FitResult.chi2min,TwinqUqUfrac.FitResult.chi2min,Twinall.FitResult.chi2min];%TwinallbutqUfrac.FitResult.chi2min,
    
    dof = Twin.FitResult.dof;
    save(save_file,'mNuSq','mNuSqErr','E0','E0Err','N','NErr','B','BErr','dof','chi2min');
 end
%% plot
if strcmp(SavePlot,'ON')
Twin.GetPlotColor;
twinLabel = {'twin',sprintf('same \\rhod'),'same T_2','same qUf','same qU','same qU qUf','all same'};%,'all same -qUf'
lineArg = {'o-','LineWidth',4,'MarkerSize',10,'MarkerFaceColor',Twin.PlotColor,'Color',rgb('IndianRed')};
x = 1:numel(twinLabel);

% neutrino mass shift
fig1 = figure(1);
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
plot(x,mNuSq,lineArg{:});
xticks(1:numel(x));
xticklabels(twinLabel);
xtickangle(35)
ylabel(sprintf('m^2_\\nu shift (eV^2)'));
PrettyFigureFormat;
grid on;

% % endpoint
fig2 = figure(2);
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.7]);
plot(x,E0,lineArg{:});
xticks(1:numel(x));
xticklabels(twinLabel)
xtickangle(35)
ylabel(sprintf('E0 shift (eV)'));
PrettyFigureFormat;
grid on;

% chi2 /dof
fig3 = figure(3);
set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.6]);
plot(x,chi2min./dof,lineArg{:});
xticks(1:numel(x));
xticklabels(twinLabel)
xtickangle(35)
ylabel(sprintf('\\chi2 / %.0f dof',dof));
PrettyFigureFormat;
xlim([0.5,numel(x)+0.5])
grid on;

% save plot

    savedirplots = strrep(savedir,'results','plots');
    savedirplots = strrep(savedirplots,'knm1Twins','knm1_Stacking');
    if ~exist(savedirplots,'dir')
        system(['mkdir ',savedirplots]);
    end
    figname = [savedirplots,sprintf('TwinStackingEffect_%s_%.0feVrange',RunList,FitRange)];
    
    print(fig1,[figname,'mNu.png'],'-dpng','-r450');
    print(fig1,[figname,'E0.png'],'-dpng','-r450');
    print(fig1,[figname,'chi2.png'],'-dpng','-r450');
end
end
