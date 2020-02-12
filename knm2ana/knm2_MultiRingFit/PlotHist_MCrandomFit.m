function PlotHist_MCrandomFit(varargin)
% function to plot simple hisogram of randomized MC data
% perfrom Gaussian fit
p=inputParser;
p.addParameter('FitResults','',@(x)iscell(x)); % cell containing fit results (structs, as assigned in RunAnalysis)
p.addParameter('FitDist','Normal',@(x)ismember(x,{'OFF','Normal','Poisson'}));
p.addParameter('PlotPar','mNu',@(x)ismember(x,{'mNu','E0','Norm','Bkg','qU','mTSq'})); %parameter which is plotted
p.addParameter('Ring',2,@(x)isfloat(x)); % for multiring parameters only
p.addParameter('SaveAs','',@(x)ischar(x)); % save plot: path + name
p.parse(varargin{:})

FitResults = p.Results.FitResults;
FitDist    = p.Results.FitDist;
PlotPar    = p.Results.PlotPar;
Ring       = p.Results.Ring;
SaveAs     = p.Results.SaveAs;

% if get rid of fits that failed
FailIndex = cell2mat(cellfun(@(x) isempty(x),FitResults,'UniformOutput',false))';
FitResults = FitResults(~FailIndex);

mypar  = cell2mat(cellfun(@(x) x.par,FitResults,'UniformOutput',false))';
myerr  = cell2mat(cellfun(@(x) x.err,FitResults,'UniformOutput',false))';
mychi2 = cell2mat(cellfun(@(x) x.chi2min,FitResults,'UniformOutput',false))';
mydof  = cell2mat(cellfun(@(x) x.dof,FitResults,'UniformOutput',false))';

nPixels = (size(mypar,1)-9)/4;
switch PlotPar
    case 'mNu'
        PlotIndex = 1;
        xstr = sprintf('{\\itm}_\\nu^2 (eV^2)');
        legstr = sprintf('\\langle{\\itm}_\\nu^2\\rangle');
        legUnit = sprintf('eV^2');
        Ring = 1; %common to all rings
    case 'E0'
        PlotIndex = 2;
        xstr = sprintf('{\\itE}_0^{fit} (eV)');
        legstr = sprintf('\\langle{\\itE}_0^{fit}\\rangle');
        legUnit = sprintf('eV');
        Ring = 1; %common to all rings
        mypar(2,:) = mypar(2,:)+18573.7; % init value
    case 'Bkg'
        BkgIndex = 3:2+nPixels;
        PlotIndex = BkgIndex(Ring);
        xstr = sprintf('B ring %.0f (mcps)',Ring);
        legstr = sprintf('\\langleB\\rangle');
        legUnit = sprintf('cps');
         mypar(PlotIndex,:) = mypar(PlotIndex,:)*1e3; % init value
    case 'Norm'
        NormIndex = 3+nPixels:3+2*nPixels-1;
        PlotIndex = NormIndex(Ring);
        xstr = sprintf('N ring %.0f',Ring);
        legstr = sprintf('\\langleN\\rangle');
        legUnit = sprintf('');
    case 'qU'
        qUIndex = 2*nPixels+9:3*nPixels+8;
        PlotIndex = qUIndex(Ring);
        xstr = sprintf('qU offset ring %.0f (eV)',Ring);
        legstr = sprintf('\\langleqU offset\\rangle');
        legUnit = sprintf('eV');
    case 'mTSq'
        mTSqIndex = 3*nPixels+10:4*nPixels+9;
        PlotIndex = mTSqIndex(Ring);
        xstr = sprintf('Energy smearing \\sigma^2 %.0f (eV^2)',Ring);
        legstr = sprintf('\\langle\\sigma^2\\rangle');
        legUnit = sprintf('eV^2');
end

if ~strcmp(FitDist,'OFF')
% Fit to distribution
dFit = fitdist(mypar(PlotIndex,:)',FitDist);
end

% plot histogram
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
h1 = histfit(mypar(PlotIndex,:));
h1(1).FaceAlpha = 1;
h1(1).FaceColor = rgb('DodgerBlue');
h1(2).Color = rgb('Black'); h1(2).LineWidth = 3;
PrettyFigureFormat('FontSize',24);
xlabel(xstr);
ylabel('Events');

if ~strcmp(FitDist,'OFF')
    if strcmp(FitDist,'Normal')
        FitDist = 'Gaussian'; %nicer in leg
    end
    leg = legend(h1(2),sprintf('%s fit: %s = %.2f %s',FitDist,legstr,dFit.mu,legUnit));
    leg.Location = 'northwest';
    leg.FontSize = get(gca,'FontSize');
    legend boxoff
end

ylim([0,max(h1(1).YData)*1.2]);

if ~isempty(SaveAs)
    savename = sprintf('%s_%s_Ring%.0f.pdf',SaveAs,PlotPar,Ring);
    export_fig(f1,savename,'-painters');
    fprintf('Save plot to %s \n',savename);
end
end