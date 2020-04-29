function out = GetSmartKsn1Grid(varargin)
p = inputParser;
p.addParameter('SanityPlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('range',90,@(x)isfloat(x));
p.addParameter('nGridSteps',50,@(x)isfloat(x));
p.parse(varargin{:});
SanityPlot = p.Results.SanityPlot;
range      = p.Results.range;
nGridSteps = p.Results.nGridSteps;

savedir  = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
filename = sprintf('%sSterileTestPar_KNM1_Twin_E0BkgNorm_%.0feVrange_50GridSteps.mat',savedir,range);

if ~exist(filename,'file')
    out = 0;
    return
end

mnu4Sq              = logspace(-1,log10(range^2),nGridSteps)';

d = importdata(filename);
ApproxSin2T4Contour = interp1(d.mnu4Sq_contour,d.sin2T4_contour,mnu4Sq,'spline');
ApproxSin2T4Contour(ApproxSin2T4Contour>0.5) = 0.5;
sin2Low = ApproxSin2T4Contour-0.1; 
sin2Low(sin2Low<=1e-03) = 1e-03;

sin2Up = ApproxSin2T4Contour+0.1;
sin2Up(sin2Up>0.5) = 0.5;

sin2LowLog = (log10(sin2Low));
sin2UpLog  = (log10(sin2Up));

%sin2T4             = cell2mat(arrayfun(@(a,b) linspace(a,b,nGridSteps),sin2Low,sin2Up,'UniformOutput',0));
sin2T4             = cell2mat(arrayfun(@(a,b) logspace(a,b,nGridSteps),sin2LowLog,sin2UpLog,'UniformOutput',0));


mnu4Sq  = repmat(mnu4Sq,1,nGridSteps);


if strcmp(SanityPlot,'ON')
    GetFigure;
    plot(sin2T4,mnu4Sq(:,1),'.');
    hold on;
    plot(ApproxSin2T4Contour,mnu4Sq,'k-','LineWidth',2);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    PrettyFigureFormat;
    xlabel('|U_{e4}|^2');
    ylabel(sprintf('{\\itm}_4 (eV^2)'));
    
end


end


