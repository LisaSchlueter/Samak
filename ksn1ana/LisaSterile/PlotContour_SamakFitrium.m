% plot script for KSN1 sterile contour plot
% Lisa, April 2020
%% settings
%Nathan = 'OFF'; % compare with nathan
Fitrium = 'ON';
% Giunti = 'OFF';
% Troitsk = 'OFF';
CL = [95];%, 95 99];
SavePlot = 'ON';
PlotSplines = 'ON'; % plot smooth spline instead of points -> fails for closed contours
range = 95;%
nGridSteps = 25;
chi2Str = {'chi2Stat','chi2CMShape'};
DataType = 'Twin';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
SmartGrid = 'OFF';
SysBudget = 24;
%% load grid (or calculate if doesn't exist)

mnu4Sq   = cell(numel(chi2Str),1);
sin2T4   = cell(numel(chi2Str),1);
chi2     = cell(numel(chi2Str),1);
chi2_ref = cell(numel(chi2Str),1);

for i=1:numel(chi2Str)
    GridArg = {'range',range,...
    'nGridSteps',nGridSteps,...
    'chi2',chi2Str{i},...
    'freePar',freePar,...
    'RunList',RunList,...
    'SmartGrid',SmartGrid,...
    'RecomputeFlag','OFF',...
     'SysBudget',SysBudget};
[mnu4Sq{i},sin2T4{i},chi2{i},chi2_ref{i},savefile] = KSN1GridSearch(...
    'DataType',DataType,GridArg{:});
end
%% collect plot argument

if numel(chi2Str)==1
    if strcmp(chi2Str,'chi2Stat')
        chi2Label = 'stat. only';
    else
        chi2Label = 'stat. and syst.';
    end
    titleStr = sprintf('%s , %s , %.0f eV at %.0f%% C.L.',DataType,chi2Label,range,CL);
else
    titleStr = sprintf('%s , %.0f eV range , %.0f%% C.L.',DataType,range,CL);  
end

for i=1:numel(chi2Str)
PlotArg ={'mnu4Sq',mnu4Sq{i},...
    'sin2T4',sin2T4{i},...
    'chi2',chi2{i},'chi2_ref',chi2_ref{i},...
    'CL',CL,...
    'titleStr',titleStr,'LineStyle','-'};
%% plot
if i>1
    [p2,~] = KSN1ContourPlot(PlotArg{:},'PlotSplines',PlotSplines,'HoldOn','ON','Color','Orange');
else
    [p1,~] = KSN1ContourPlot(PlotArg{:},'PlotSplines',PlotSplines,'Color','DodgerBlue');
end

end

%% nathan
% if strcmp(Nathan,'ON')
%     if range==95
%         nrange = 90;
%     else
%         nrange = range;
%     end
%     nathanfile = sprintf('%sNathan/NathanContour_%.0feV_%s_%s.mat',savedir,nrange,DataType,chi2Str);
%     legStr = [legStr,{'Nathan'}];
%     
%     try
%         d = importdata(nathanfile);
%         hold on;
%         pN = plot(d.sith4_X,d.m4_Y,'-.','LineWidth',2);
%     catch
%         fprintf('nathan file %s not available \n',nathanfile)
%     end
% end

savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
    %% fitrum
    if strcmp(Fitrium,'ON')

      for i=1:numel(chi2Str)
        fitriumfile = sprintf('%sOthers/contour_KSN1_Fitrium_%s_%.0feV_%s_%.0f_0.txt',...
            savedir,DataType,range,chi2Str{i},CL);
        try
            d = importdata(fitriumfile);
            hold on;
            if i==1
                pFstat = plot(d.data(:,1),d.data(:,2),'-.','LineWidth',2,'Color',p1.Color);
            elseif i==2
                pFsys = plot(d.data(:,1),d.data(:,2),'-.','LineWidth',2,'Color',p2.Color);
            end
        catch
            fprintf('fitrium file %s not available \n',fitriumfile)
        end
      end
     
    end
%% legend

if numel(chi2Str)==1
    legStr = {'Samak','Fitrium'};
    leg = legend(legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
else
    legStr = {'Samak (stat. only)','Fitrium (stat. only)','Samak (stat. and syst.)','Fitrium (stat. and syst.)'};
    leg = legend([p1,pFstat,p2,pFsys],legStr,'EdgeColor',rgb('Silver'),'Location','southwest');
end


leg.Title.String = 'Analysis';
leg.Title.FontWeight = 'normal';
 %xlim([3e-03 0.5])

%% save plot

%% plot label
plotdir = [getenv('SamakPath'),'ksn1ana/LisaSterile/plots/'];
MakeDir(plotdir);
plotname = sprintf('%sKSN1_CompareFitrium_%.0feV_%s.png',plotdir,range,DataType);

if strcmp(SavePlot,'ON')
    print(plotname,'-dpng','-r450');
    fprintf('save plot to %s \n',plotname);
end


%% 
% giunti
% if strcmp(Giunti,'ON')
%     legStr = [legStr,{'Giunti'}];
%     giuntifile = sprintf('%sOthers/coord_Giunti.mat',savedir);
%     try
%         dG = importdata(giuntifile);
%         hold on;
%         pGiunti = plot(dG.sith4_X.^2,dG.m4_Y.^2,'-.','LineWidth',2);
%     catch
%         fprintf('giunti file %s not available \n',giuntifile)
%     end
% end
% 
% % troitsk
% if strcmp(Troitsk,'ON')
%     legStr = [legStr,{'Troitsk'}];
%     troitskfile = sprintf('%sOthers/coord_troitsk.mat',savedir);
%     try
%         dT = importdata(troitskfile);
%         hold on;
%         pT = plot(dT.sith4_X.^2,dT.m4_Y.^2,'--','LineWidth',2);
%     catch
%         fprintf('trotsk file %s not available \n',nathanfile)
%     end
%     
% end