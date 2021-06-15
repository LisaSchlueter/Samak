% plot alternative pixel list results
nFits = 500;
freePar = 'mNu E0 Bkg Norm';
DataType = 'Real';%'Real';
Parameter = 'mNuSq';
RunList = 'KNM2_Prompt';
range = 40;                % fit range in eV below endpoint
AltPixLists = {'Half','AziHalfEW','AziHalfNS'};
savedir = [getenv('SamakPath'),'knm2ana/knm2_AltRunPixLists/results/'];

%% load random half
savenameRand = sprintf('%sknm2_PixListRandHalf_%s_%s_%.0feV_%.0ffits_KNM2.mat',...
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
nPixel = zeros(numel(AltPixLists),2);

for i=1:numel(AltPixLists)
    savename = sprintf('%sknm2_PixListAlt_%s_%s_%s_%.0feV_KNM2.mat',...
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
        nPixel(i,:) = cellfun(@(x) numel(x),d.PixList)';
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
f1 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.5]);
meanY =mean(y);%wmean(y,1./yErr.^2);
h1 =histogram(y-meanY,'FaceColor',rgb('LightGray'),'FaceAlpha',1,'Normalization','probability');
xlabel(sprintf('%s - \\langle%s\\rangle (%s)',xStr,xStr,Unit));
ylabel('Frequency');
PrettyFigureFormat('FontSize',22)
titleStr = sprintf('\\langle%s\\rangle = %.2f %s , \\sigma =  %.3f %s',xStr,meanY,Unit,std(y),Unit);
t = title(sprintf('%.0f stacked runs - %s',numel(RunList),titleStr),...
    'FontWeight','normal','FontSize',get(gca,'FontSize')-2);
legStrRand = sprintf('Random %.0f pixels',size(PixList,2));
hold on;
yLimits = ylim;
%% plot alternative 
yplot = linspace(0,1000,20);
legStr = {legStrRand};
Sigma = cell(numel(AltPixLists),1);

for i=1:numel(AltPixLists)
     AltY = ParAlt{i};
  
    
    Sigma{i} = (AltY-meanY)./std(y);
    
    if strcmp(AltPixLists{i},'Half')
        legStr = [legStr,sprintf('Inner, %.0f pixel: %s - \\langle%s\\rangle   = %.2f eV^{ 2}',nPixel(i,1),xStr,xStr,AltY(1)-meanY),...
           sprintf('Outer, %.0f pixel: %s - \\langle%s\\rangle  = %.2f eV^{ 2}',nPixel(i,2),xStr,xStr,AltY(2)-meanY)];
        PlotArg = {'LineWidth',3,'Color',rgb('Orange')};
        LineStyle = {'-',':'};
    elseif strcmp(AltPixLists{i},'Azi')
        legStr = [legStr,'Azi 1','Azi 2','Azi 3','Azi 4'];
        PlotArg = {'LineWidth',3,'Color',rgb('MediumSeaGreen')};
        LineStyle = {'-','-.','--',':'};
    elseif strcmp(AltPixLists{i},'AziHalfEW')
        strcmp(AltPixLists{i},'Half')
        legStr = [legStr,...
            sprintf('East, %.0f pixel: %s - \\langle%s\\rangle    = %.2f eV^{ 2}',nPixel(i,1),xStr,xStr,AltY(1)-meanY),...
            sprintf('West, %.0f pixel: %s - \\langle%s\\rangle   = %.2f eV^{ 2}',nPixel(i,2),xStr,xStr,AltY(2)-meanY)];
        PlotArg = {'LineWidth',3,'Color',rgb('MediumSeaGreen')};
        LineStyle = {'-',':'};
    elseif strcmp(AltPixLists{i},'AziHalfNS')
        strcmp(AltPixLists{i},'Half')
        legStr = [legStr,...
            sprintf('North, %.0f pixel: %s - \\langle%s\\rangle  = %.2f eV^{ 2}',nPixel(i,1),xStr,xStr,AltY(1)-meanY),...
            sprintf('South, %.0f pixel: %s - \\langle%s\\rangle  = %.2f eV^{ 2}',nPixel(i,2),xStr,xStr,AltY(2)-meanY)];
        PlotArg = {'LineWidth',3,'Color',rgb('DodgerBlue')};
        LineStyle = {'-',':'};
    end
    
     for n=1:nSeg(i)
        plot(AltY(n).*ones(20,1)-meanY,yplot,LineStyle{n},PlotArg{:})
    end
end

leg = legend(legStr{:},'Location','northwest');
leg.EdgeColor = rgb('Silver');
leg.NumColumns = 2;
ylim([0,yLimits(2)+0.1]);
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.9])); 

xmin = (min(cell2mat(ParAlt))-meanY);
xmax = (max(cell2mat(ParAlt))-meanY);
xlim(1.5*[xmin,xmax])

plotname = strrep(strrep(savenameRand,'results','plots'),'.mat','_PlotRandAlt.pdf');
export_fig(f1,plotname);
fprintf('save plot to %s \n',plotname)
%% standard deviations from mean
Sigma_v = cell2mat(Sigma);
