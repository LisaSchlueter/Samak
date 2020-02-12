% function to make nice correlation plot with histograms
function CorrHistPlot(par1,par2,FitResults,varargin)
%%
p=inputParser;
p.addParameter('SaveAs','',@(x)ischar(x)); % save plot: path + name
p.addParameter('Ring',2,@(x)isfloat(x)); % for multiring parameters only
p.addParameter('PlotCorMat','ON',@(x)ismember(x,{'ON','OFF'})); % for multiring parameters only

p.parse(varargin{:})
SaveAs     = p.Results.SaveAs;
Ring       = p.Results.Ring;
PlotCorMat = p.Results.PlotCorMat;

nPixels = (size(FitResults,1)-9)/4;
if strcmp(par1,'mNu')
    x = FitResults(1,:);
    xstr = sprintf('{\\itm_\\nu}^2 (eV^{2})');
elseif strcmp(par1,'E0')
    x = FitResults(2,:);
    xstr =sprintf('{\\itE}^{fit}_0 (eV)');
elseif strcmp(par1,'N')
    x = FitResults(4+nPixels,:)+1; %  normalization
    xstr =sprintf('{\\itN}_{%.0f} ',Ring);
elseif strcmp(par1,'qU')
    x = FitResults(2*nPixels+10,:); % qU-offset
    xstr =sprintf('\\Delta{\\itqU}_{%.0f} (eV)',Ring);
elseif strcmp(par1,'mTSq')
    x = FitResults(3*nPixels+11,:); % mTSq
    xstr =sprintf('\\sigma_{%.0f}^{2} (eV^2)',Ring);
end

if strcmp(par2,'mNu')
    y = FitResults(1,:);
    ystr = sprintf('{\\itm_\\nu}^2 (eV^{2})');
elseif strcmp(par2,'E0')
    y = FitResults(2,:);
    ystr =sprintf('{\\itE}^{fit}_0 (eV)');
elseif strcmp(par2,'N')
    y = FitResults(4+nPixels,:)+1; % normalization
    ystr =sprintf('{\\itN}_{%.0f}',Ring);
elseif strcmp(par2,'qU')
    y = FitResults(2*nPixels+10,:); %  qU-offset
    ystr =sprintf('\\Delta{\\itqU}_{%.0f} (eV)',Ring);
    strcmp(par1,'mTSq')
elseif strcmp(par2,'mTSq')
    y = FitResults(3*nPixels+11,:); %  mTSq
    ystr =sprintf('\\sigma^{%.0f} (eV^2)',Ring);
end

xcentral = mean(x);
ycentral = mean(y);

f111 = figure('Renderer','opengl');
set(f111, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
p = scatterhist(x,y,'Direction','out',...);%,'HistogramDisplayStyle','bar',...
    'Location','Northeast','Color',rgb('DodgerBlue'));
hold on;
xlabel(xstr);
ylabel(ystr);
PrettyFigureFormat('FontSize',28);
thisYlim = ylim; thisXlim =xlim;
plot(xcentral.*[1,1],[ycentral,max(y)*5],':','Color',rgb('Orange'),'LineWidth',3);
plot([xcentral,max(x)*5],ycentral.*[1,1],':','Color',rgb('Orange'),'LineWidth',3);
pfit = plot(xcentral,ycentral,'d','MarkerSize',12,'Color',rgb('DarkGoldenRod'),...
    'MarkerFaceColor',rgb('Orange'),'LineWidth',2);
xlim(thisXlim); ylim(thisYlim);
leg = legend(sprintf('%.0f Pseudo experiments',numel(x)),'Mean');
leg.Location='northwest';
legend boxoff;

if ~isempty(SaveAs)
    savename = sprintf('%s_%s%sCorr_Ring%.0f_ScatterPlot.pdf',SaveAs,par1,par2,Ring);
    export_fig(f111,savename,'-painters');
    fprintf('Save plot as %s \n',savename);
end

if strcmp(PlotCorMat,'ON')
    %% correlation matrix all fit parameters
    ParIndex = FitResults(:,1)~=0;
    FitResultsActive = FitResults(ParIndex,:);
    if nPixels<=2
        f112 = figure('Units','normalized','Position',[0.1,0.1,0.6,0.6]);
    elseif nPixels==4
        f112 = figure('Units','normalized','Position',[0.1,0.1,0.9,1]);
    end
    CorrMat = corrcoef(FitResultsActive');
    cp= corplot(CorrMat);
    if nPixels==2
        myticklabels = {sprintf('{\\itm}_\\nu^2'),sprintf('{\\itE}_{0}^{fit}'),...
            sprintf('{\\itB_1}'),sprintf('{\\itB_2}'),sprintf('{\\itN_1}'),sprintf('{\\itN_2}'),sprintf('\\DeltaqU_2'),sprintf('\\sigma_2^2')};
    elseif nPixels==4
        myticklabels = {sprintf('{\\itm}_\\nu^2'),sprintf('{\\itE}_{0}^{fit}'),...
            sprintf('{\\itB_1}'),sprintf('{\\itB_2}'),sprintf('{\\itB_3}'),sprintf('{\\itB_4}')...
            sprintf('{\\itN_1}'),sprintf('{\\itN_2}'),sprintf('{\\itN_3}'),sprintf('{\\itN_4}'),...
            sprintf('\\DeltaqU_2'),sprintf('\\DeltaqU_3'),sprintf('\\DeltaqU_4'),...
            sprintf('\\sigma_2^2'),sprintf('\\sigma_3^2'),sprintf('\\sigma_4^2')};
    end
    xticklabels(myticklabels);
    yticklabels(myticklabels);
    
    cp.LineWidth= get(gca,'LineWidth');
    c = colorbar;
    c.Label.String = sprintf('correlation coefficient \\rho');
    c.Label.FontSize = get(gca,'FontSize');
    c.LineWidth = get(gca,'LineWidth');
    
    if nPixels <=2
        PrettyFigureFormat('FontSize',25);
        c.FontSize = get(gca,'FontSize')-2;
        toffset = 0.15;
    elseif nPixels==4
        PrettyFigureFormat('FontSize',20);
        c.FontSize = get(gca,'FontSize')+4;
        toffset = 0.2;
    end
    
    for i=1:size(CorrMat,1)
        for j=1:size(CorrMat,1)
            text(i+toffset,j+0.5,sprintf('%.2f',CorrMat(i,j)),'FontSize',get(gca,'FontSize')-3);
        end
    end
    set(gca,'XMinorTick','off');
    set(gca,'YMinorTick','off');
    set(gca,'TickLength',[0 0]);
    
    
    if ~isempty(SaveAs)
        savename = sprintf('%s_%s%sCorr_Ring%.0f_Matrix.pdf',SaveAs,par1,par2,Ring);
        export_fig(f112,savename,'-painters');
        fprintf('Save plot as %s \n',savename);
    end
end
end
