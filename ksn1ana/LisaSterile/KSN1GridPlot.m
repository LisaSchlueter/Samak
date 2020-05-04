function KSN1GridPlot(varargin)
p = inputParser;
p.addParameter('mnu4Sq','',@(x)isfloat(x));
p.addParameter('sin2T4','',@(x)isfloat(x));
p.addParameter('chi2','',@(x)isfloat(x));
p.addParameter('chi2_ref','',@(x)isfloat(x));
p.addParameter('CL',0.9,@(x)isfloat(x));
p.addParameter('SaveAs','OFF',@(x) ischar(x));
p.addParameter('titleStr','',@(x)ischar(x));
p.addParameter('nInter',1e3,@(x)isfloat(x));
p.addParameter('ContourPlot','OFF',@(x)ismember(x,{'ON','OFF','Fitrium'}));
p.parse(varargin{:});

mnu4Sq      = p.Results.mnu4Sq;
sin2T4      = p.Results.sin2T4;
chi2        = p.Results.chi2;
chi2_ref    = p.Results.chi2_ref;
CL          = p.Results.CL;
SaveAs      = p.Results.SaveAs;
titleStr    = p.Results.titleStr;
nInter      = p.Results.nInter;
ContourPlot = p.Results.ContourPlot;

if isempty(mnu4Sq)
    load('KSN1_SmartGridInit_40eVrange.mat','mnu4Sq','sin2T4','chi2','chi2_ref');
    fprintf('no input grid specified - load default grid \n')
    titleStr = 'Default grid 40 eV range (twins)';
end
%% 

%% mesh grid interpolation -> even finer binning, nicer appearance!
nGridSteps = size(mnu4Sq,1);

[X,Y] = meshgrid(mnu4Sq(:,1),sin2T4(1,:));
mNugrid    = repmat(logspace(log10(min(min(mnu4Sq))),log10(max(max(mnu4Sq))),nInter),nInter,1);
sin2T4grid = repmat(logspace(log10(min(min(sin2T4))),log10(max(max(sin2T4))),nInter),nInter,1)';
chi2grid = reshape(interp2(X,Y,chi2,mNugrid,sin2T4grid),nInter ,nInter );

% chi2 = chi2grid';
% sin2T4 = sin2T4grid;
% mnu4Sq = mNugrid;
%%
GetFigure;
DeltaChi2 = GetDeltaChi2(CL,2);
chi2_ref = min(min(chi2grid));

chi2grid((chi2grid-chi2_ref)>DeltaChi2) =  NaN;%DeltaChi2+chi2_ref;% NaN;
zlimMax = DeltaChi2;
%zlimMax = DeltaChi2;
surf(sin2T4grid,mNugrid,chi2grid-chi2_ref,'EdgeColor','interp','FaceColor','interp');
zlim([0 zlimMax])
set(gca,'XScale','log')
set(gca,'YScale','log')
c =colorbar;
c.Label.String = sprintf('\\Delta\\chi^2');
PrettyFigureFormat('FontSize',24)
c.Label.FontSize = get(gca,'FontSize')+2;
%c.Limits=[0 zlimMax];
xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
zlabel(sprintf('\\Delta\\chi^2'))
view([0 0 1])
grid off


if strcmp(ContourPlot,'ON')
    [mnu4Sq_contour, sin2T4_contour] = ...
        KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,CL);
    y = linspace(min(mnu4Sq_contour),max(mnu4Sq_contour)+1e3,1e3);
    x = interp1(mnu4Sq_contour,sin2T4_contour,y,'spline');
    hold on;
    plot3(x,y,DeltaChi2.*ones(numel(x),1),'k-','LineWidth',2);
    SaveAs = strrep(SaveAs,'.pdf','_C.pdf');
    SaveAs = strrep(SaveAs,'.png','_C.png');
    % plot3(sin2T4_contour,mnu4Sq_contour,zeros(numel(sin2T4_contour),1));
elseif strcmp(ContourPlot,'Fitrium')
    fitriumfile0 = 'contour_KSN1_Fitrium_Real_95eV_chi2CMShape_95_0.txt';
    fitriumfile1 = 'contour_KSN1_Fitrium_Real_95eV_chi2CMShape_95_1.txt';
    fitriumfile2 = 'contour_KSN1_Fitrium_Real_95eV_chi2CMShape_95_2.txt';
    d0 = importdata(fitriumfile0);
    d1 = importdata(fitriumfile1);
    d2 = importdata(fitriumfile2);
    hold on; 
   plot3(d0.data(:,1),d0.data(:,2),DeltaChi2.*ones(numel(d0.data(:,1)),1),'k-','LineWidth',2);
   hold on;
   plot3(d1.data(:,1),d1.data(:,2),DeltaChi2.*ones(numel(d1.data(:,1)),1),'k-','LineWidth',2);
   plot3(d2.data(:,1),d2.data(:,2),DeltaChi2.*ones(numel(d2.data(:,1)),1),'k-','LineWidth',2);
   SaveAs = strrep(SaveAs,'.png','_FitriumContour.png');
     SaveAs = strrep(SaveAs,'.pdf','_FitriumContour.pdf');
end

rangeApprox = sqrt(max(max(mnu4Sq)));
if rangeApprox<90
    xlim([5e-03 0.5])
elseif rangeApprox>=90
    xlim([1e-03 0.5])
end
ylim([1 (rangeApprox)^2]);
title(titleStr,'FontWeight','normal','FontSize',get(gca,'FontSize'))

if ~strcmp(SaveAs,'OFF')
    PlotDir = [getenv('SamakPath'),'ksn1ana/LisaSterile/plots/'];
    MakeDir(PlotDir);
    SaveAs = [PlotDir,SaveAs];
    if contains(SaveAs,'.pdf')
        export_fig(gcf,SaveAs);
    elseif contains(SaveAs,'.png')
        print(gcf,SaveAs,'-dpng','-r450');
    end
    fprintf('save plot to %s \n',SaveAs)
end
end