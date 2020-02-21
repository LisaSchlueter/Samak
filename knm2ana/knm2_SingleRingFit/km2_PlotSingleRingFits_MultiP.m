% plot Single Ring fits for all RW periods

%% settings
RunList = [1,2,3];
nRuns   = [121,95,92];
freePar = 'mNu E0 Bkg Norm';
Range   = 90;
ROIFlag =  '14keV';%Default';%'14keV';%'Default';
RingMerge = 'Full';
chi2 = 'chi2Stat';

PlotPar      = 1;
SavePlot     = 'ON';
linFitFlag   = 'ON';
PlotMode     = 'Abs';
Blind        = 'ON';
YLim         = '';%[-0.28,0.18];[-11 5];%%

if PlotPar==1
    PlotMode = 'Rel'; % Blinding
end
%% retreive / calculate results
switch ROIFlag
    case 'Default'
        RoiStr = '';
    case '14keV'
        RoiStr = '_14keVROI';
end

savedir = [getenv('SamakPath'),'tritium-data/fit/Knm2/SingleRingFit/'];
savename = arrayfun(@(x,y) sprintf('%sSingleRingFitResult_%s_KNM2_RW%.0f_%.0fruns_fix%s_%s_%.0feVrange%s.mat',...
    savedir,RingMerge,x,y,strrep(freePar,' ',''),...
    chi2,Range,RoiStr),RunList,nRuns,...
    'UniformOutput',false);

if all(cellfun(@(x) exist(x,'file'),savename))
    d = cellfun(@(x) importdata(x),savename,'UniformOutput',false);
    fprintf('Load results from file %s \n',savename{1})
else
    d = cell(numel(savename),1);
    CommonArg = {'freePar',freePar,'Range',Range,'ROIFlag',ROIFlag,'chi2',chi2,'RingMerge',RingMerge};
    InputArg = {CommonArg{:},'RunList','KNM2_RW1';...
        CommonArg{:},'RunList','KNM2_RW2';...
        CommonArg{:},'RunList','KNM2_RW3'};
    
    parfor i=1:numel(d)
        if exist(savefile{i},'file')
            d{i} = importdata(savefile{i});
        else
            knm2_SingleRingFit(InputArg{i,:});
            d{i} = importdata(savefile{i});
        end
    end
end

%% plot

AllPar = cell2mat(cellfun(@(x) x.par,d,'UniformOutput',false)');
AllPar(:,4) = AllPar(:,4)+1;% normalization
AllErr = cell2mat(cellfun(@(x) x.err,d,'UniformOutput',false)');
AllErrNeg = cell2mat(cellfun(@(x) x.errNeg,d,'UniformOutput',false)');
AllErrPos = cell2mat(cellfun(@(x) x.errPos,d,'UniformOutput',false)');
nPeriods = numel(d);
nRings = size(AllPar,1)/nPeriods;
RingList = 1:nRings;

y       = zeros(nPeriods,nRings);
yErr    = zeros(nPeriods,nRings);
yErrNeg = zeros(nPeriods,nRings);
yErrPos = zeros(nPeriods,nRings);

for i=1:nPeriods
    StartI    = ((i-1)*nRings)+1;
    EndI      = i*nRings;
    y(i,:)    = AllPar(StartI:EndI,PlotPar);
    yErr(i,:) = AllErr(StartI:EndI,PlotPar);
    yErrNeg(i,:) = AllErrNeg(StartI:EndI,PlotPar);
    yErrPos(i,:) = AllErrPos(StartI:EndI,PlotPar);
end

if PlotPar==1 % neutrino mass: use average MINOS errors
   meanErr = (yErrPos - yErrNeg)./2;
   yErr(meanErr~=0) = meanErr(meanErr~=0);
end
%% linear fit
linFitpar = zeros(nPeriods,2);
linFiterr = zeros(nPeriods,2);
linFitchi2min = zeros(nPeriods,1);
linFitdof = zeros(nPeriods,1);
if strcmp(linFitFlag,'ON')
    for i=1:nPeriods
        [linFitpar(i,:), linFiterr(i,:), linFitchi2min(i),linFitdof(i)] = linFit(RingList',y(i,:)',yErr(i,:)');
    end
end

%% plot fit result
fig2 = figure('Renderer','painters');
set(fig2,'units','normalized','pos',[0.1, 0.1,0.6,0.5]);
plot(linspace(0.5,nRings+0.5,10),zeros(10,1),'-','Color',rgb('SlateGrey'),'LineWidth',2);
hold on;

if strcmp(PlotMode,'Rel')
    meanPar =wmean(y',1./yErr'.^2)';
elseif strcmp(PlotMode,'Abs')
    meanPar = zeros(nPeriods,1);
end
Colors = {'DodgerBlue','GoldenRod','IndianRed'};
Ebar = cell(nPeriods,1);
l    = cell(nPeriods,1);
EbarArg = {'o','LineWidth',2,'LineStyle','none','MarkerSize',8};
LineStyles = {'-',':','-.'};
for i=1:nPeriods
    Ebar{i} = errorbar(RingList,(y(i,:)-meanPar(i)),yErr(i,:),EbarArg{:},...
        'Color',rgb(Colors{i}),'MarkerFaceColor',rgb(Colors{i}));
    Ebar{i}.CapSize = 0;
    if strcmp(linFitFlag,'ON')
       l{i} = plot(RingList,(linFitpar(i,1).*RingList+linFitpar(i,2)-meanPar(i))',...
           'Color',rgb(Colors{i}),'LineWidth',2.5,'LineStyle',LineStyles{i});
    end
end
PrettyFigureFormat('FontSize',24);
%% legend
switch PlotPar
    case 1
        UnitStr = sprintf('eV^2');
        Prefix = 'm';
    case 2
        UnitStr = sprintf('eV');
          Prefix = 'm';
    case 4
        UnitStr = sprintf('');
        Prefix = sprintf('10^{-3}');
end
if strcmp(linFitFlag,'ON')
    if any(abs(linFitpar(:,1))<0.1)
        linFitleg1 =  sprintf('RW 1: (%.1f \\pm %.1f) %s%s/ring , p-value = %.2f',linFitpar(1,1)*1e3,linFiterr(1,1)*1e3,Prefix,UnitStr,1-chi2cdf(linFitchi2min(1),linFitdof(1)));
        linFitleg2 =  sprintf('RW 2: (%.1f \\pm %.1f) %s%s/ring , p-value = %.2f',linFitpar(2,1)*1e3,linFiterr(2,1)*1e3,Prefix,UnitStr,1-chi2cdf(linFitchi2min(2),linFitdof(2)));
        linFitleg3 =  sprintf('RW 3: (%.1f \\pm %.1f) %s%s/ring , p-value = %.2f',linFitpar(3,1)*1e3,linFiterr(3,1)*1e3,Prefix,UnitStr,1-chi2cdf(linFitchi2min(3),linFitdof(3)));
    else
        linFitleg1 =  sprintf('RW 1: (%.1f \\pm %.1f) %s/ring , p-value = %.2f',linFitpar(1,1),linFiterr(1,1),UnitStr,1-chi2cdf(linFitchi2min(1),linFitdof(1)));
        linFitleg2 =  sprintf('RW 2: (%.1f \\pm %.1f) %s/ring , p-value = %.2f',linFitpar(2,1),linFiterr(2,1),UnitStr,1-chi2cdf(linFitchi2min(2),linFitdof(2)));
        linFitleg3 =  sprintf('RW 3: (%.1f \\pm %.1f) %s/ring , p-value = %.2f',linFitpar(3,1),linFiterr(3,1),UnitStr,1-chi2cdf(linFitchi2min(3),linFitdof(3)));
    end
    leg = legend([l{:}],linFitleg1,linFitleg2,linFitleg3);
    %sprintf('%s , %.0f eV range',runname,Range),
else
    leg = legend('RW 1','RW 2','RW 3');
end
leg.EdgeColor = rgb('Silver');
%leg.Title.String = sprintf('%.0f eV range',Range);
hold off;
leg.FontSize= get(gca,'FontSize')-2;
leg.Location = 'best';
%% axis labels
switch PlotMode
    case 'Rel'
        if PlotPar==1
            ylabel(sprintf('{\\itm}_\\nu^2 - \\langle{\\itm}_\\nu^2\\rangle (eV^2)'));
        elseif PlotPar==2
            ylabel(sprintf('{\\itE}_0^{fit} - \\langle{\\itE}_0^{fit}\\rangle (eV)'));
        elseif  PlotPar==4
            ylabel(sprintf('{\\itN} - \\langle{\\itN}\\rangle '));
        end
    case 'Abs'
        if PlotPar==1
            ylabel(sprintf('{\\itm}_\\nu^2 (eV^2)'));
        elseif PlotPar==2
            ylabel(sprintf('{\\itE}_0^{fit} - 18573.7 (eV)'));
        elseif  PlotPar==4
            ylabel(sprintf('{\\itN} '));
        end
end
xlabel('Rings')

switch RingMerge
    case 'InnerPseudo3'
        xticks([1]);
        xticklabels({'1,2,3,4,5,6,7,8,9'})
    case 'Full'
        xticks([1 2 3 4]);
        xticklabels({'1-3','4-6','7-9','10-12'})
    case 'Default'
        xticks([1:1:10]);
        xticklabels({'1-2',3,4,5,6,7,8,9,10,'11-12'})
    case 'None'
        xticks([1:12])
end
set(gca,'XMinorTick','off');
xlim([min(RingList)-0.2,max(RingList)+0.2])

ymin = min(min((y-meanPar)-yErr));
ymax = max(max((y-meanPar)+yErr));
if ~isempty(YLim)
    ylim([min(YLim),max(YLim)]);
elseif PlotPar==4
    ylim([ymin-5e-03,ymax+5e-03])
elseif strcmp(PlotMode,'Rel') || ymin<0
    ylim([1.1*ymin,1.9*ymax]);
else
    ylim([0.5*ymin,1.5*ymax]);
end

%% title
t = title(sprintf('%.0f eV range, %s ROI, Fit parameter: %s ',Range,ROIFlag,freePar));
t.FontWeight = 'normal';
t.FontSize = get(gca,'FontSize');

%% save
if strcmp(SavePlot,'ON')
    savedir = [getenv('SamakPath'),'tritium-data/plots/Knm2/SingleRingFit/'];
    MakeDir(savedir);
    plotname = strrep(strrep(strrep(savename{1},'.mat',sprintf('_%.0f.pdf',PlotPar)),'KNM2_RW1_121runs_',''),'fit/','plots/');
    export_fig(gcf,plotname);
    fprintf('Save plot to %s \n',plotname);
end