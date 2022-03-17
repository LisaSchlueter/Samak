
% plot alternative and random run lists
% KNM2, stacked uniform 40 eV range
AltRunLists = {'KNM2_Up','KNM2_Down','KNM2_RW1','KNM2_RW2','KNM2_RW3'};       % defines alternative pixel list
nAltRunLists = numel(AltRunLists);
freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';
range = 40;               % fit range in eV below endpoint
Parameter = 1;            % 1 = numass, 2= endpoint
Mode = 'Rel';
FSDFlag = 'KNM2_0p1eV';
BKG_PtSlope = 3*1e-06;
%% load random half run lists fit results
nFits = 2;
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];
savenameRand = sprintf('%sknm2_RunListRandHalf_%s_%s_%.0feV_%.0ffits_%s_BkgPT%.2gmucpsPers.mat',...
    savedir,DataType,strrep(freePar,' ',''),range,nFits,FSDFlag,BKG_PtSlope*1e6);
dRand = importdata(savenameRand);
ParRand = cellfun(@(x) x.par(Parameter),dRand.FitResults);
ErrRand = cellfun(@(x) x.err(Parameter),dRand.FitResults);
switch Parameter
    case 1
        ErrNegRand = cellfun(@(x) x.errNeg(Parameter),dRand.FitResults);
        ErrPosRand = cellfun(@(x) x.errPos(Parameter),dRand.FitResults);
        ErrRand    = 0.5*(ErrPosRand-ErrNegRand);
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
ParAlt   = zeros(nAltRunLists,1);
ErrAlt   = zeros(nAltRunLists,1);
nRunsAlt = zeros(nAltRunLists,1);
FilesAlt = cell(nAltRunLists,1);
for i=1:nAltRunLists
   savenameAlt = sprintf('%sknm2_AltRunList_%s_%s_%s_%.0feV_%s_BkgPt%.2g.mat',...
    savedir,AltRunLists{i},DataType,strrep(freePar,' ',''),range,FSDFlag,BKG_PtSlope*1e6);
    FilesAlt{i} = importdata(savenameAlt);
    switch Parameter
        case 1
            ParAlt(i) =  FilesAlt{i}.mNuSq;
            ErrAlt(i) =  0.5*(FilesAlt{i}.FitResult.errPos(1)-FilesAlt{i}.FitResult.errNeg(1));
        case 2
            ParAlt(i) =  FilesAlt{i}.E0;
            ErrAlt(i) =  FilesAlt{i}.E0Err;
    end
    nRunsAlt(i) = numel( FilesAlt{i}.RunList);
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
legStrRand = sprintf('Random %.0f runs (%.0f samples)',numel(dRand.RunLists{1}),nFits);
hold on;
yLimits = ylim;

%% plot alternative 
yplot = linspace(0,1000,20);
legStr = {legStrRand};
Sigma = zeros(nAltRunLists,1);

for i=1:numel(AltRunLists)
    
    if strcmp(AltRunLists{i},'KNM2_Up')
        legStr = [legStr,sprintf('Up scans (%.0f runs):      {\\itm}_\\nu^2 = (%.2f \\pm %.2f) eV^2',nRunsAlt(i),ParAlt(i),ErrAlt(i))];
        PlotArg = {'LineWidth',3,'Color',rgb('Orange')};
        LineStyle = {'-'};
    elseif strcmp(AltRunLists{i},'KNM2_Down')
        legStr = [legStr,sprintf('Down scans (%.0f runs):  {\\itm}_\\nu^2 = (%.2f \\pm %.2f) = %.2f eV^2',nRunsAlt(i),ParAlt(i),ErrAlt(i))];
        PlotArg = {'LineWidth',3,'Color',rgb('Crimson')};
        LineStyle = {':'};
    elseif contains(AltRunLists{i},'KNM2_RW')
        RWperiod = str2num(extractAfter(AltRunLists{i},'RW'));
        if RWperiod==1
            Urw = -49.6; % mV
                  PlotArg = {'LineWidth',3,'Color',rgb('DodgerBlue')};
        LineStyle = {'--'};
        elseif RWperiod==2
             Urw = -7.7; % mV
                   PlotArg = {'LineWidth',3,'Color',rgb('LimeGreen')};
        LineStyle = {'-.'};
        else
             Urw = 193; % mV
                   PlotArg = {'LineWidth',3,'Color',rgb('Navy')};
        LineStyle = {'-'};
        end
             legStr = [legStr,sprintf('P%.0f: {\\itU}_{RW} = %.0f mV (%.0f runs):  (%.2f \\pm %.2f)  = %.2f eV^2',...
                 RWperiod,Urw,nRunsAlt(i),ParAlt(i),ErrAlt(i))];    
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

