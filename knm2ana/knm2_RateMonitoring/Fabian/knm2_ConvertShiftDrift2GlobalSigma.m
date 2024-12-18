function SigmaGlobal = knm2_ConvertShiftDrift2GlobalSigma(varargin)
% convert individual shifts + drifts/broadenings per rear wall period into 1 global sigma
% can be applied to uniform or pseudo-ring analysis
p = inputParser;
p.addParameter('sigmas_v',[0.05,0.03,0.08],@(x)isfloat(x));  % broadenings for 3 periods
p.addParameter('shifts_v',[0.2 -0.1 0],@(x)isfloat(x));      % shifts for 3 periods
p.addParameter('weights_v',[171,93,97]./361,@(x)isfloat(x)); % relative time spent in each period
p.addParameter('SanityPlot','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('PlotName','./plots/Knm2_GlobalDrift.pdf',@(x)ischar(x));
p.parse(varargin{:});

sigmas_v   = p.Results.sigmas_v;
shifts_v   = p.Results.shifts_v;
weights_v  = p.Results.weights_v;
SanityPlot = p.Results.SanityPlot;
PlotName   = p.Results.PlotName;

sigmas_v = sigmas_v.*1e3;
shifts_v = shifts_v.*1e3;
% grid to define pdf
Emin = min(shifts_v)-5.*max(sigmas_v);
Emax = max(shifts_v)+5.*max(sigmas_v);
E = linspace(Emin,Emax,1e3);

GaussFun       = @(x,mu,sigma,weight) weight./(sigma*sqrt(2*pi))*exp(-0.5*(x-mu).^2./sigma.^2);
GlobalGaussFun = @(x) GaussFun(x,shifts_v(1),sigmas_v(1),weights_v(1))+...
                   GaussFun(x,shifts_v(2),sigmas_v(2),weights_v(2))+...
                   GaussFun(x,shifts_v(3),sigmas_v(3),weights_v(3));
               
% binned pdf and cdf               
%GlobalGauss    = GlobalGaussFun(E);    
GlobalGaussCDF = cumsum(GlobalGaussFun(E))./max(cumsum(GlobalGaussFun(E)));
[GlobalGaussCDF, Index] = unique(GlobalGaussCDF,'stable');
E = E(Index);

% get pdf samples from inverse cdf
nSamples = 1e5;
GlobalGaussSamples = interp1(GlobalGaussCDF,E,rand(nSamples,1),'lin');
GlobalGaussSamples(abs(GlobalGaussSamples)>1e3)= [];

SigmaGlobal = std(GlobalGaussSamples);

% sanity plot
if strcmp(SanityPlot,'ON')
    figure('Units','normalized','Position',[0.1,0.1,0.5,0.4])
    h1=histogram(GlobalGaussSamples,'FaceColor',rgb('DodgerBlue'),'Normalization','pdf');
    hold on;
    plot(E,GlobalGaussFun(E),'LineWidth',2,'Color',rgb('Orange'));%.*h1.BinWidth.*nSamples
    hold off  
    PrettyFigureFormat('FontSize',22)
    leg = legend(sprintf('Global \\sigma_t = %.0f meV',SigmaGlobal),'Effective source potential model: 3 Gaussians');
    PrettyLegendFormat(leg);
    xlim([min(shifts_v)-2*max(sigmas_v),max(shifts_v)+2*max(sigmas_v)]);
    ylim([0 0.01]);
    xlabel('Source potential (meV)');
    ylabel('Probability density');
 %   print(gcf,PlotName,'-dpng','-r500');
    export_fig(PlotName);
end
