function pd1 = FitEOdistribution(arg_1)
%CREATEFIT    Create plot of datasets and fits
%   PD1 = CREATEFIT(ARG_1)
%   Creates a plot, similar to the plot in the main distribution fitter
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1
%
%   See also FITDIST.

% This function was automatically generated on 24-May-2018 09:19:58

% Output fitted probablility distribution: PD1

% Data from dataset "results(:,2) data":
%    Y = arg_1 (originally results(:,2))

% Force all inputs to be column vectors
arg_1 = arg_1(:);

% Prepare figure
clf;
hold on;
LegHandles = []; LegText = {};


% --- Plot data originally in dataset "results(:,2) data"
[CdfF,CdfX] = ecdf(arg_1,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 3;
BinInfo.nbins = 10;
[~,BinEdge] = internal.stats.histbins(arg_1,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
    'LineStyle','-', 'LineWidth',1);
xlabel('Effective Endpoint -18575 (eV)');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'Toy MC simulations';

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);


% --- Create fit "fit 1"

% Fit this distribution to get parameter values
% To use parameter estimates from the original fit:
%     pd1 = ProbDistUnivParam('normal',[ -0.5411564518306, 2.408231484418])
pd1 = fitdist(arg_1, 'normal');
YPlot = pdf(pd1,XGrid);
hLine = plot(XGrid,YPlot,'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
LegHandles(end+1) = hLine;
LegText{end+1} = 'Gaussian Fit';

% Adjust figure
box on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
set(hLegend,'Interpreter','none');
