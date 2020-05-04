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
p.parse(varargin{:});

mnu4Sq      = p.Results.mnu4Sq;
sin2T4      = p.Results.sin2T4;
chi2        = p.Results.chi2;
chi2_ref    = p.Results.chi2_ref;
CL          = p.Results.CL;
SaveAs      = p.Results.SaveAs;
titleStr    = p.Results.titleStr;
nInter      = p.Results.nInter;
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
%chi2_ref = chi2(1,1); 
chi2grid((chi2grid-chi2_ref)>DeltaChi2) =  NaN;%DeltaChi2+chi2_ref;% NaN;
%chi2((chi2-chi2_ref)<0) =  NaN;% -1+chi2_ref;% NaN;
surf(sin2T4grid,mNugrid,chi2grid-chi2_ref,'EdgeColor','interp','FaceColor','interp');
zlim([0 DeltaChi2])
set(gca,'XScale','log')
set(gca,'YScale','log')
c =colorbar;
c.Label.String = sprintf('\\Delta\\chi^2');
PrettyFigureFormat('FontSize',24)
c.Label.FontSize = get(gca,'FontSize')+2;
c.Limits=[0 DeltaChi2];
xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
zlabel(sprintf('\\Delta\\chi^2'))
view([0 0 1])
 grid off
 
 ContourPlot = 'OFF';
if strcmp(ContourPlot,'ON')
     [mnu4Sq_contour, sin2T4_contour] = ...
        KSN1Grid2Contour(mnu4Sq,sin2T4,chi2,chi2_ref,CL);
    hold on;
    surf(mnu4Sq_contour, sin2T4_contour,zeros(numel(sin2T4_contour),1));
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