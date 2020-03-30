% plot alternative pixel list results
nFits = 500;
freePar = 'mNu E0 Bkg Norm';
DataType = 'Twin';%'Real';
Parameter = 'E0';
RunList = 'KNM2_Prompt';
range = 40;                % fit range in eV below endpoint
AltPixLists = {'Half','AziHalfEW','AziHalfNS'};
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];

%% load random half
savenameRand = sprintf('%sknm2_PixListRandHalf_%s_%s_%.0feV_%.0ffits.mat',...
    savedir,DataType,strrep(freePar,' ',''),range,nFits);
if exist(savenameRand,'file')
    load(savenameRand);
else
    fprintf(2,'Random half pix list not computed yet \n')
    return
end
%% load Alt run lists
ParAlt = cell(numel(AltPixLists),1);
nSeg   = zeros(numel(AltPixLists),1);

for i=1:numel(AltPixLists)
    savename = sprintf('%sknm2_PixListAlt_%s_%s_%s_%.0feV.mat',...
        savedir,AltPixLists{i},DataType,strrep(freePar,' ',''),range);
    
    if exist(savename,'file')
        d = importdata(savename);
        switch Parameter
            case 'mNuSq'
                ParAlt{i} = d.mNuSq;
            case 'E0'
                ParAlt{i} = d.E0;
        end
        nSeg(i) = numel(d.E0);
    else
        fprintf('alternative pix list %s not computed yet \n',AltPixLists{i})
    end
end
%% get parameters
switch Parameter
    case 'mNuSq'
        y    = mNuSq;
        yErr = mNuSqErr;
        xStr = sprintf('{\\itm}_\\nu^2');
        Unit = sprintf('eV^2');
    case 'E0'
        y    = E0+Q_i;
        yErr = E0Err;
        xStr = sprintf('{\\itE}_0^{fit}');
        Unit = sprintf('eV');
end
%% plot random half
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
meanY =wmean(y,1./yErr.^2);
h1 =histogram(y-meanY,'FaceColor',rgb('LightGray'),'FaceAlpha',1);
xlabel(sprintf('%s - \\langle%s\\rangle (%s)',xStr,xStr,Unit));
ylabel('Occurrence');
PrettyFigureFormat('FontSize',22)
titleStr = sprintf('\\langle%s\\rangle = %.2f %s , \\sigma =  %.3f %s',xStr,meanY,Unit,std(y),Unit);
t = title(sprintf('%.0f stacked runs - %s',numel(RunList),titleStr),...
    'FontWeight','normal');
legStrRand = sprintf('Random %.0f pixels',size(PixList,2));
hold on;
yLimits = ylim;
%% plot alternative 
yplot = linspace(0,1000,20);
legStr = {legStrRand};
Sigma = cell(numel(AltPixLists),1);

for i=1:numel(AltPixLists)
    
    if strcmp(AltPixLists{i},'Half')
        legStr = [legStr,'Inner','Outer'];
        PlotArg = {'LineWidth',3,'Color',rgb('Orange')};
        LineStyle = {'-',':'};
    elseif strcmp(AltPixLists{i},'Azi')
        legStr = [legStr,'Azi 1','Azi 2','Azi 3','Azi 4'];
        PlotArg = {'LineWidth',3,'Color',rgb('MediumSeaGreen')};
        LineStyle = {'-','-.','--',':'};
    elseif strcmp(AltPixLists{i},'AziHalfEW')
        strcmp(AltPixLists{i},'Half')
        legStr = [legStr,'East','West'];
        PlotArg = {'LineWidth',3,'Color',rgb('MediumSeaGreen')};
        LineStyle = {'-',':'};
    elseif strcmp(AltPixLists{i},'AziHalfNS')
        strcmp(AltPixLists{i},'Half')
        legStr = [legStr,'North','South'];
        PlotArg = {'LineWidth',3,'Color',rgb('DodgerBlue')};
        LineStyle = {'-',':'};
    end
    
    AltY = ParAlt{i};
    for n=1:nSeg(i)
        plot(AltY(n).*ones(20,1)-meanY,yplot,LineStyle{n},PlotArg{:})
    end
    
    Sigma{i} = (AltY-meanY)./std(y);
end

leg = legend(legStr{:},'Location','northwest');
leg.EdgeColor = rgb('Silver');
ylim(yLimits);

xmin = (min(cell2mat(ParAlt))-meanY);
xmax = (max(cell2mat(ParAlt))-meanY);
xlim(1.5.*[xmin,xmax])

plotname = strrep(strrep(savenameRand,'results','plots'),'.mat','_PlotRandAlt.pdf');
export_fig(f1,plotname);
fprintf('save plot to %s \n',plotname)
%% standard deviations from mean
Sigma_v = cell2mat(Sigma);
