function KSN1GridPlot(varargin)
p = inputParser;
p.addParameter('mnu4Sq','',@(x)isfloat(x));
p.addParameter('sin2T4','',@(x)isfloat(x));
p.addParameter('chi2','',@(x)isfloat(x));
p.addParameter('chi2_ref','',@(x)isfloat(x));
p.addParameter('CL',0.9,@(x)isfloat(x));
p.addParameter('SaveAs','OFF',@(x) ischar(x));
p.addParameter('titleStr','',@(x)ischar(x));
p.parse(varargin{:});

mnu4Sq      = p.Results.mnu4Sq;
sin2T4      = p.Results.sin2T4;
chi2        = p.Results.chi2;
chi2_ref    = p.Results.chi2_ref;
CL          = p.Results.CL;
SaveAs      = p.Results.SaveAs;
titleStr    = p.Results.titleStr;

if isempty(mnu4Sq)
    load('KSN1_SmartGridInit_40eVrange.mat','mnu4Sq','sin2T4','chi2','chi2_ref');
    fprintf('no input grid specified - load default grid \n')
    titleStr = 'Default grid 40 eV range (twins)';
end

% PlotContour = 'OFF';
% range =95;%
% nGridSteps =50;%25;%50;
% chi2Str = 'chi2CMShape';
% DataType = 'Real';
% freePar = 'E0 Bkg Norm';
% RunList = 'KNM1';
% SmartGrid = 'OFF';
% CL = 0.82;
% [mnu4Sq,sin2T4,chi2,chi2_ref,savefile] = KSN1GridSearch('range',range,...
%     'nGridSteps',nGridSteps,...
%     'chi2',chi2Str,...
%     'DataType',DataType,...
%     'freePar',freePar,...
%     'RunList',RunList,...
%     'SmartGrid',SmartGrid,...
%     'RecomputeFlag','OFF');

%
%%
GetFigure;
DeltaChi2 = GetDeltaChi2(CL,2);
chi2_ref = min(min(chi2));
%chi2_ref = chi2(1,1); 
chi2((chi2-chi2_ref)>DeltaChi2) =  DeltaChi2+chi2_ref;% NaN;
%chi2((chi2-chi2_ref)<0) =  NaN;% -1+chi2_ref;% NaN;
surf(sin2T4,mnu4Sq,chi2'-chi2_ref,'EdgeColor','interp','FaceColor','interp');
zlim([0 DeltaChi2])
set(gca,'XScale','log')
set(gca,'YScale','log')
c =colorbar;
c.Label.String = sprintf('\\Delta\\chi^2');
PrettyFigureFormat('FontSize',24)
c.Label.FontSize = get(gca,'FontSize')+2;
c.Limits=[0 DeltaChi2];
xlabel('|U_{e4}|^2');
ylabel(sprintf('{\\itm}_4 (eV^2)'));
zlabel(sprintf('\\Delta\\chi^2'))
view([0 0 1])
 grid off
 
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