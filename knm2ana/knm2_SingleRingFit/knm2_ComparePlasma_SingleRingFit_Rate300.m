% plot Single Ring fits for all RW periods

%% settings
RunList = [1,2,3];
nRuns   = [121,95,92];
freePar = 'E0 Bkg Norm';
Range   = 90;
ROIFlag =  '14keV';%Default'
RingMerge = 'Full';
chi2 = 'chi2Stat';
MosCorrFlag = 'OFF';
SavePlot     = 'ON';
YLim         = [-150,100];%'';%[-0.28,0.18];[-11 5];%%

%% retreive / calculate results
switch ROIFlag
    case 'Default'
        RoiStr = '';
    case '14keV'
        RoiStr = '_14keVROI';
end
if strcmp(MosCorrFlag,'ON')
    MosStr = '_MosCorr';
else
    MosStr = '';
end
savedir = [getenv('SamakPath'),'tritium-data/fit/Knm2/SingleRingFit/'];
savename = arrayfun(@(x,y) sprintf('%sSingleRingFitResult_%s_KNM2_RW%.0f_%.0fruns_fix%s_%s_%.0feVrange%s%s.mat',...
    savedir,RingMerge,x,y,strrep(freePar,' ',''),...
    chi2,Range,RoiStr,MosStr),RunList,nRuns,...
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

%% prepare variables
PlotPar = 2;  %endpoint
AllPar = cell2mat(cellfun(@(x) x.par,d,'UniformOutput',false)');
AllPar(:,4) = AllPar(:,4)+1;% normalization
AllErr = cell2mat(cellfun(@(x) x.err,d,'UniformOutput',false)');
AllErrNeg = cell2mat(cellfun(@(x) x.errNeg,d,'UniformOutput',false)');
AllErrPos = cell2mat(cellfun(@(x) x.errPos,d,'UniformOutput',false)');
nPeriods = numel(d);
nRings = size(AllPar,1)/nPeriods;
RingList = 1:nRings;

E0       = zeros(nPeriods,nRings);
E0Err    = zeros(nPeriods,nRings);

for i=1:nPeriods
    StartI    = ((i-1)*nRings)+1;
    EndI      = i*nRings;
    E0(i,:)    = AllPar(StartI:EndI,PlotPar);
    E0Err(i,:) = AllErr(StartI:EndI,PlotPar);
end

%% load 300ev analysis results
if strcmp(MosCorrFlag,'ON')
    yRM =  [64,72,-78;...
           46,23,-112;...
          32,-3,-134;...
         -8,--58,-179]'.*1e-3;
    
    % 300eV analysis is relative -> set to reference Period2, ring1 reference
    yRM = yRM+(E0(2,1)-yRM(2,1));
else
    yRM =  [86,40,-104;...
            67,-9,-138;...
            52,-34,-160;...
            13,-90,-203]'.*1e-3;
end

y = (E0-yRM)*1e3;
%y2 = yRM*1e3;
yErr = E0Err*1e3;
%% plot fit result
fig2 = figure('Renderer','painters');
set(fig2,'units','normalized','pos',[0.1, 0.1,0.6,0.5]);
plot(linspace(0.5,nRings+0.5,10),zeros(10,1),'-','Color',rgb('SlateGrey'),'LineWidth',2);
hold on;
Colors = {'DodgerBlue','GoldenRod','IndianRed'};
Ebar = cell(nPeriods,1);
pRM  = cell(nPeriods,1);
l    = cell(nPeriods,1);
EbarArg = {'o','LineWidth',2,'LineStyle','none','MarkerSize',8};
LineStyles = {'-',':','-.'};
for i=1:nPeriods
    Ebar{i} = errorbar(RingList,y(i,:),yErr(i,:),EbarArg{:},...
        'Color',rgb(Colors{i}),'MarkerFaceColor',rgb(Colors{i}));
    Ebar{i}.CapSize = 0;
%     pRM = plot(RingList+0.1,yRM(i,:)*1e3,'d','LineWidth',2,'LineStyle','none','MarkerSize',8,...
%          'Color',rgb(Colors{i}),'MarkerFaceColor',rgb(Colors{i}));
end
PrettyFigureFormat('FontSize',24);
%% legend

leg = legend([Ebar{:}],'RW 1','RW 2','RW 3');
leg.EdgeColor = rgb('Silver');
%leg.Title.String = sprintf('%.0f eV range',Range);
hold off;
leg.FontSize= get(gca,'FontSize')-2;
leg.Location = 'best';
%% axis labels
ylabel(sprintf('{\\itE}_0^{fit} - U_{300 eV} (meV)'))
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

ymin = min(min(y)-E0Err);
ymax = max(max(y)+E0Err);
if ~isempty(YLim)
    ylim([min(YLim),max(YLim)]);
% elseif ymin<0
%     ylim([1.1*ymin,1.9*ymax]);
% else
%     ylim([0.5*ymin,1.5*ymax]);
end

%% title
if strcmp(MosCorrFlag,'ON')
    MosStr = ', MoS drift corrected';
else
    MosStr = '';
end
t = title(sprintf('%.0f eV range, %s ROI, Fit parameter: %s %s',Range,ROIFlag,freePar,MosStr));
t.FontWeight = 'normal';
t.FontSize = get(gca,'FontSize');

%% save
if strcmp(SavePlot,'ON')
    savedir = [getenv('SamakPath'),'tritium-data/plots/Knm2/SingleRingFit/'];
    MakeDir(savedir);
    plotname = strrep(strrep(strrep(savename{1},'.mat',sprintf('_%.0f_ComparePlasmaShift300eV.pdf',PlotPar)),'KNM2_RW1_121runs_',''),'fit/','plots/');
    export_fig(gcf,plotname);
    fprintf('Save plot to %s \n',plotname);
end