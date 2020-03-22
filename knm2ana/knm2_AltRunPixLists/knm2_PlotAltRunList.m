
% plot alternative and random run lists
% KNM2, stacked uniform 40 eV range
AltRunLists = {'KNM2_Up','KNM2_Down'};       % defines alternative pixel list
nAltRunLists = numel(AltRunLists);
freePar = 'E0 Bkg Norm';
DataType = 'Real';
range = 40;               % fit range in eV below endpoint
Parameter = 2;            % 1 = numass, 2= endpoint
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

meanParRand = wmean(ParRand,1./ErrRand.^2);
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
h1 =histogram(ParRand-meanParRand,'FaceColor',rgb('LightGray'),'FaceAlpha',1);
xlabel(sprintf('%s - \\langle%s\\rangle (%s)',xStr,xStr,Unit));
ylabel('Occurrence');
PrettyFigureFormat('FontSize',22)
titleStr = sprintf('\\langle%s\\rangle = %.2f %s , \\sigma =  %.3f %s',...
    xStr,meanParRand,Unit,std(ParRand),Unit);
t = title(sprintf('Uniform fit - %s',titleStr),...
    'FontWeight','normal');
legStrRand = sprintf('Random %.0f runs',numel(dRand.RunList{1}));
hold on;
yLimits = ylim;

%% plot alternative 
yplot = linspace(0,1000,20);
legStr = {legStrRand};
Sigma = zeros(nAltRunLists,1);

for i=1:numel(AltRunLists)
    
    if strcmp(AltRunLists{i},'KNM2_Up')
        legStr = [legStr,sprintf('Up scans (%.0f runs)',nRunsAlt(i))];
        PlotArg = {'LineWidth',3,'Color',rgb('Orange')};
        LineStyle = {'-'};
    elseif strcmp(AltRunLists{i},'KNM2_Down')
        legStr = [legStr,sprintf('Down scans (%.0f runs)',nRunsAlt(i))];
        PlotArg = {'LineWidth',3,'Color',rgb('Crimson')};
        LineStyle = {':'};
    end
    
     plot(ParAlt(i).*ones(20,1)-meanParRand,yplot,LineStyle{:},PlotArg{:})
    Sigma(i) = (ParAlt(i)-meanParRand)./std(ParRand);
end

leg = legend(legStr{:},'Location','northwest');
leg.EdgeColor = rgb('Silver');
ylim(yLimits);

plotname = strrep(strrep(savenameRand,'results','plots'),'.mat','_PlotRandAlt.pdf');
export_fig(f1,plotname);
fprintf('save plot to %s \n',plotname);
%xmin = (min(cell2mat(ParAlt))-meanParRand);
%xmax = (max(cell2mat(ParAlt))-meanParRand);
%xlim(1.5.*[xmin,xmax])

