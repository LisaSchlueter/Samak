
% plot alternative and random run lists
% KNM2, stacked uniform 40 eV range
AltRunLists = {'KNM2_Up','KNM2_Down'};       % defines alternative pixel list
nAltRunLists = numel(AltRunLists);
freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;               % fit range in eV below endpoint
Parameter = 1;            % 1 = numass, 2= endpoint
Mode = 'Rel';
%% load random half run lists fit results
nFits = 500;
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savenameRand = sprintf('%sknm2_RunListRandHalf_%s_%s_%.0feV_%.0ffits.mat',...
    savedir,DataType,strrep(freePar,' ',''),range,nFits);
dRand = importdata(savenameRand);
ParRand = cellfun(@(x) x.par(Parameter),dRand.FitResult);
ErrRand = cellfun(@(x) x.err(Parameter),dRand.FitResult);
switch Parameter
    case 1
        xStr = sprintf('{\\itm}_\\nu^2');
        Unit = sprintf('eV^2');
    case 2
        Q_i = 18573.7;
        ParRand = ParRand+Q_i;
        xStr = sprintf('{\\itE}_0^{fit}');
        Unit = sprintf('eV');
end

% if strcmp(Mode,'Abs')
%     meanParRand =0;
% else
meanParRand = mean(ParRand);%wmean(ParRand,1./ErrRand.^2);
%end
%% load alternative runs lists
ParAlt  = zeros(nAltRunLists,1);
nRunsAlt = zeros(nAltRunLists,1);
for i=1:nAltRunLists
    savenameAlt = sprintf('%sknm2_AltRunList_%s_%s_%s_%.0feV.mat',...
        savedir,AltRunLists{i},DataType,strrep(freePar,' ',''),range);
    d = importdata(savenameAlt);
    switch Parameter
        case 1
            ParAlt(i) = d.mNuSq;
        case 2
            ParAlt(i) = d.E0;
    end
    nRunsAlt(i) = numel(d.RunList);
end

%% plot random half
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
h1 =histogram(ParRand-meanParRand,'FaceColor',rgb('LightGray'),'FaceAlpha',1,'Normalization','probability');
xlabel(sprintf('%s - \\langle%s\\rangle (%s)',xStr,xStr,Unit));
ylabel('Frequency');
PrettyFigureFormat('FontSize',20)
titleStr = sprintf('\\langle%s\\rangle = %.2f %s , \\sigma =  %.2f %s',...
    xStr,meanParRand,Unit,std(ParRand),Unit);
t = title(sprintf('%s',titleStr),...
    'FontWeight','normal','FontSize',get(gca,'FontSize')-1);
legStrRand = sprintf('Random %.0f runs (%.0f samples)',numel(dRand.RunList{1}),nFits);
hold on;
yLimits = ylim;

%% plot alternative 
yplot = linspace(0,1000,20);
legStr = {legStrRand};
Sigma = zeros(nAltRunLists,1);

for i=1:numel(AltRunLists)
    
    if strcmp(AltRunLists{i},'KNM2_Up')
        legStr = [legStr,sprintf('Up scans, %.0f runs:      {\\itm}_\\nu^2 - \\langle{\\itm}_\\nu^2\\rangle = %.2f eV^2',nRunsAlt(i),ParAlt(i))];
        PlotArg = {'LineWidth',3,'Color',rgb('Orange')};
        LineStyle = {'-'};
    elseif strcmp(AltRunLists{i},'KNM2_Down')
        legStr = [legStr,sprintf('Down scans, %.0f runs:  {\\itm}_\\nu^2 - \\langle{\\itm}_\\nu^2\\rangle = %.2f eV^2',nRunsAlt(i),ParAlt(i))];
        PlotArg = {'LineWidth',3,'Color',rgb('Crimson')};
        LineStyle = {':'};
    end
    
     plot(ParAlt(i).*ones(20,1)-meanParRand,yplot,LineStyle{:},PlotArg{:})
    Sigma(i) = (ParAlt(i)-meanParRand)./std(ParRand);
end

leg = legend(legStr{:},'Location','northwest');
leg.EdgeColor = rgb('Silver');
leg.FontSize = get(gca,'FontSize');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.85])); 
ylim([yLimits(1),yLimits(2)+0.06]);

plotname = strrep(strrep(savenameRand,'results','plots'),'.mat','_PlotRandAlt.pdf');
export_fig(f1,plotname);
fprintf('save plot to %s \n',plotname);
%xmin = (min(cell2mat(ParAlt))-meanParRand);
%xmax = (max(cell2mat(ParAlt))-meanParRand);
%xlim(1.5.*[xmin,xmax])

