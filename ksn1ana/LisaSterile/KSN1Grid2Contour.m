function [mnu4sq_contour, sin2t4_contour] = KSN1Grid2Contour(mnu4sq,sin2t4,chi2,chi2_ref,CL)
%% calculate contour from msq4 - sin2t4 grid at confidence level cl

% define delta chi2 for given confidence level
DeltaChi2 = GetDeltaChi2(CL,2);

nGridStep = size(mnu4sq,1);
sin2t4_contour = zeros(nGridStep,1);
mnu4sq_contour = mnu4sq(:,1);

if chi2_ref>min(min(chi2))
     chi2_ref = min(min(chi2));
end
for i=1:nGridStep
    if max(chi2(:,i)-chi2_ref)<DeltaChi2 || min(chi2(:,i)-chi2_ref)>DeltaChi2
        sin2t4_contour(i) = NaN;
%    else %chi2(1,i)>min(chi2(:,i))
%         % points with local minium    
%        IndexMin = find(chi2(:,i)==min(chi2(:,i)),1);    % find minimum
%        %[sinInter, IndexInter] = unique(sin2t4(i,IndexMin:end));
%         [chi2Inter, IndexInter] = unique(chi2(IndexMin:end,i));
%         sinInter = sin2t4(IndexMin:end,i);
%        sin2t4_contour(i) = interp1(chi2Inter-chi2_ref,sinInter(IndexInter),DeltaChi2,'spline');
     else
         sin2t4_contour(i) = interp1(chi2(:,i)-chi2_ref,sin2t4(i,:),DeltaChi2,'spline');
    end
end

ExclLogic = ( abs(sin2t4_contour)>=1) | isnan(sin2t4_contour);
sin2t4_contour(ExclLogic) = [];
mnu4sq_contour(ExclLogic) = [];

nGridStep = numel(mnu4sq_contour);
%% sanity check
SanityPlot = 'OFF';%OFF';
if strcmp(SanityPlot,'ON')
    GetFigure;
    mnuIndex = 15;
    tStr1 = 'problem m4 =';
    tStr2 = 'good m4 =';
    for i=1:nGridStep
        if any((chi2(:,i)-chi2_ref)<-0.1)
            PlotArg = {'-.','Color',rgb('FireBrick'),'LineWidth',2};
            tStr1 = [tStr1,sprintf(' %.0f ,',sqrt(mnu4sq_contour(i)))];
        else
            PlotArg = {'-','Color',rgb('DimGray'),'LineWidth',2};
             tStr2 = [tStr2,sprintf(' %.0f ,',sqrt(mnu4sq_contour(i)))];
        end
        p = plot(sin2t4(i,:),chi2(:,i)-chi2_ref, PlotArg{:});
        xlabel(sprintf('|U_{e4}|^2'));
        ylabel(sprintf('\\Delta \\chi^2'));
        PrettyFigureFormat;
        hold on;
    end
    tStr = sprintf('%s eV \n%s eV',tStr1,tStr2);
    title(tStr,'FontWeight','normal','FontSize',get(gca,'FontSize'));
    ylim([-20 10])
    xlim([0 0.2])
end
end