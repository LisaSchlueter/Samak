function [mnu4Sq,sin2T4] = GetSmartKsn1Grid(varargin)
p = inputParser;
p.addParameter('SanityPlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('range',95,@(x)isfloat(x));
p.addParameter('nGridSteps',100,@(x)isfloat(x));
p.addParameter('ConfLevel',90,@(x)isfloat(x));
p.addParameter('AddSin2T4',1e-02,@(x)isfloat(x)); % scan around expected contour within this region 
p.parse(varargin{:});
SanityPlot = p.Results.SanityPlot;
range      = p.Results.range;
nGridSteps = p.Results.nGridSteps;
ConfLevel  = p.Results.ConfLevel;
AddSin2T4  = p.Results.AddSin2T4;

savedir  = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/Ksn1SmartGrid/'];
filename = sprintf('%sKSN1_SmartGridInit_%.0feVrange.mat',savedir,range);

if ~exist(filename,'file')
    fprintf(2,'file does not exist %s \n',filename)
    return
end

mnu4Sq              = logspace(-1,log10((range+5)^2),nGridSteps)';

d = importdata(filename);
[mnu4Sq_contour, sin2T4_contour] = GetSterileContour(d.mnu4Sq,d.sin2T4,d.chi2,d.chi2_ref,ConfLevel);
ExclLogic = abs(sin2T4_contour)>0.5;
sin2T4_contour(ExclLogic) = [];
mnu4Sq_contour(ExclLogic)=[];

ApproxSin2T4Contour = interp1(mnu4Sq_contour,sin2T4_contour,mnu4Sq,'spline');
ApproxSin2T4Contour(ApproxSin2T4Contour>0.5) = 0.5;
sin2Low = ApproxSin2T4Contour-AddSin2T4; 
sin2Low(sin2Low<=1e-03) = 1e-03;

sin2Up = ApproxSin2T4Contour+AddSin2T4;
sin2Up(sin2Up>0.5) = 0.5;

sin2LowLog = (log10(sin2Low));
sin2UpLog  = (log10(sin2Up));

%sin2T4             = cell2mat(arrayfun(@(a,b) linspace(a,b,nGridSteps),sin2Low,sin2Up,'UniformOutput',0));
sin2T4             = cell2mat(arrayfun(@(a,b) logspace(a,b,nGridSteps),sin2LowLog,sin2UpLog,'UniformOutput',0));


mnu4Sq  = repmat(mnu4Sq,1,nGridSteps);


if strcmp(SanityPlot,'ON')
    GetFigure;
    ps = plot(sin2T4,mnu4Sq(:,1),'.');
    hold on;
    pex = plot(ApproxSin2T4Contour,mnu4Sq(:,1),'k-','LineWidth',2);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    PrettyFigureFormat;
    xlabel('|U_{e4}|^2');
    ylabel(sprintf('{\\itm}_4 (eV^2)'));
    xlim([8*1e-04, 0.7])
    legend([pex;ps(1)],'Approx. expected contour','New grid search area',...
        'EdgeColor',rgb('Silver'),'Location','southwest');
end


end


