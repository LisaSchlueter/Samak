function PlotMTD(varargin)
%PlotMTD('TDprefix','FlatSignal120','RunTimeSec',6*3600)

% --------------------------- PARSER START ------------------------------%
p = inputParser;
p.addParameter('E0',18575,@(x)isfloat(x));
p.addParameter('TDprefix','FlatSignal30');
p.addParameter('RunTimeSec',0);
p.addParameter('Save','OFF',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});
E0           = p.Results.E0;              % eV
TDprefix     = p.Results.TDprefix;        % eV
RunTimeSec   = p.Results.RunTimeSec;      % seconds
Save         = p.Results.Save;            %
% --------------------------- PARSER ENDS  ------------------------------%

% MTD Name
TD   = [TDprefix '.mat'];
MTD  = load(TD);

% Plot
figure(1)

if  RunTimeSec>0
    stairs(MTD.qU-E0, MTD.qUfrac.*RunTimeSec);
else
    stairs(MTD.qU-E0, MTD.qUfrac);
end

plt = Plot();
plt.LineWidth = .1;
ptlstr=sprintf('MTD: %s - %0.f bins',MTD.TD,numel(MTD.qU));
plt.Title = ptlstr; % plot title
MyXlabel   = sprintf('qU - %0.1f (V)',E0);
plt.XLabel = MyXlabel; % xlabel
plt.YLabel = 'Time Fraction'; %ylabel
plt.YScale = 'log';
plt.FontSize = 16;
ptlfile = sprintf('%s.png',MTD.TD);
ptl.Legend = {'MTD envelop','MTD Time Fraction'};
ptl.LegendLoc = {'NorthWest'};
hold on
%
if  RunTimeSec>0
    bar(MTD.qU-E0, MTD.qUfrac.*RunTimeSec,'blue')
else
    bar(MTD.qU-E0, MTD.qUfrac,'blue');
end

hold off
plt.YLim = [1e-5 .1];
plt.export(ptlfile);

end