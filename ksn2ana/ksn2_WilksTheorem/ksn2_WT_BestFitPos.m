% plot Asimov sensitivity + best fits on randomized contour
Hypothesis = 'H0';
SavePlt = 'OFF';
NrandMC = 1e3;
switch Hypothesis
    case 'H0'
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        chi2 = 'chi2CMShape';
    case 'H1'
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2Stat';
end
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_%.0fsamples.mat',savedir,NrandMC);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_%.0fsamples.mat',savedir,Twin_mNu4Sq,Twin_sin2T4,NrandMC);
end

if exist(savefile,'file') 
    fprintf('load file %s \n',savefile);
    d = importdata(savefile);
else
     fprintf('file missing %s \n',savefile);
    return
end

nClosed = sum(d.ClosedLog95);
IdxClosed = find(d.ClosedLog95);
PltColors = jet(nClosed);
%%
GetFigure
PltColors2 = jet(1e3);
pSens = plot(d.sin2T4_contour_Asimov,d.mNu4Sq_contour_Asimov,'-k','LineWidth',2);
hold on;
x = linspace(1e-03,0.5,1e2);
y = linspace(0.1,40^2,1e2);
if strcmp(Hypothesis,'H0')
    pref = plot(x,0.1.*ones(100,1),'-','Color',rgb('Silver'),'LineWidth',1.5);
    plot(x,40^2.*ones(100,1),'-','Color',rgb('Silver'),'LineWidth',1.5);
    plot(1e-03.*ones(100,1),y,'-','Color',rgb('Silver'),'LineWidth',1.5);
    plot(0.5.*ones(100,1),y,'-','Color',rgb('Silver'),'LineWidth',1.5);
end
% look if best fit is in our outside
InOutIdx = zeros(1e3,1);

for i=1:numel(d.chi2_bf)
     p1 = plot(d.sin2T4_bf(i),d.mNu4Sq_bf(i),'.',...
        'Color',rgb('Orange'),'MarkerSize',13,'LineWidth',1.5);

   if i~=IdxClosed  
       % delta_chi2 < 5.99
        p1.Color = rgb('SkyBlue');
        p1.LineWidth = 1;
        pHandleOpen = p1;
   else % delta_chi2 > 5.99
       pHandleReg = p1;
       p1.Color = rgb('Orange');
   end
   
   if strcmp(Hypothesis,'H0')
       %look if best fit is in our outside
       sinTmp = interp1(d.mNu4Sq_contour_Asimov,d.sin2T4_contour_Asimov,d.mNu4Sq_bf(i),'spline');
       InOutIdx(i) = sinTmp<=d.sin2T4_bf(i); % 1 means inside (significant)
       if InOutIdx(i)
           p1.Color = rgb('IndianRed');
           pHandleIn = p1;
       end
   elseif strcmp(Hypothesis,'H1')
%        [Dist, IdxDist] = sort(abs(d.mNu4Sq_bf(i)-d.mNu4Sq_contour_Asimov));
%        Diff_tmp = abs(diff(IdxDist));
%        Idx2 = find(Diff_tmp>1,1)+1;
%        LowerSin = min([d.sin2T4_contour_Asimov(IdxDist(1)),d.sin2T4_contour_Asimov(IdxDist(Idx2))]);
%        UpperSin = max([d.sin2T4_contour_Asimov(IdxDist(1)),d.sin2T4_contour_Asimov(IdxDist(Idx2))]);
%        if  d.sin2T4_bf(i) > LowerSin &&  d.sin2T4_bf(i) < UpperSin
%            
%        else
%            InOutIdx(i) = 1;
%            p1.Color = rgb('IndianRed');
%            pHandleIn = p1;
%        end  
   end
end

set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);
 if strcmp(Hypothesis,'H0')
    leg = legend([pSens,pref,pHandleOpen,pHandleReg,pHandleIn],...
        'Sensitivity','Grid boundaries',sprintf('Rand MC: \\Delta\\chi^2 < 5.99 (%.1f%%)',100.*(numel(d.sin2T4_bf)-numel(IdxClosed))./numel(d.sin2T4_bf)),...
        sprintf('Rand MC: \\Delta\\chi^2 \\geq 5.99 (%.1f%%)',100*numel(IdxClosed)./numel(d.sin2T4_bf)),...
        sprintf('Rand MC: Best fit inside sensitivity curve (%.1f%%)',100.*sum(InOutIdx)./numel(d.sin2T4_bf)),...
        'Location','southwest');
 else
      leg = legend([pSens,pHandleOpen,pHandleReg],...
        'Sensitivity',sprintf('Rand MC: \\Delta\\chi^2 < 5.99 (%.1f%%)',100.*(numel(d.sin2T4_bf)-numel(IdxClosed))./numel(d.sin2T4_bf)),...
        sprintf('Rand MC: \\Delta\\chi^2 \\geq 5.99 (%.1f%%)',100*numel(IdxClosed)./numel(d.sin2T4_bf)),...
        'Location','southwest');
%      leg = legend([pSens,pref,pHandleOpen,pHandleReg,pHandleIn],...
%         'Sensitivity','Grid boundaries',sprintf('Rand MC: \\Delta\\chi^2 < 5.99 (%.1f%%)',100.*(numel(d.sin2T4_bf)-numel(IdxClosed))./numel(d.sin2T4_bf)),...
%         sprintf('Rand MC: \\Delta\\chi^2 \\geq 5.99 (%.1f%%)',100*numel(IdxClosed)./numel(d.sin2T4_bf)),...
%         sprintf('Rand MC: Best fit outside sensitivity curve (%.1f%%)',100.*sum(InOutIdx)./numel(d.sin2T4_bf)),...
%         'Location','southwest');
end
PrettyLegendFormat(leg);
xlim([8e-04 0.65]);
ylim([0.07 50^2]);

if strcmp(SavePlt,'ON')
    plotnameContourBf = strrep(strrep(savefile,'results','plots'),'.mat','_BestFitPosition.png');
    print(gcf,plotnameContourBf,'-dpng','-r450');
    fprintf('save plot to %s \n',plotnameContourBf);
end


if strcmp(Hypothesis,'H0') && ~isfield(d,'InOutIdx')
    save(savefile,'InOutIdx','-append')
end
%%
%sum(d.sin2T4_bf==0.5)
IdxSinMax1 = find(d.sin2T4_bf==0.5 & d.mNu4Sq_bf<100);
IdxSinMax2 = find(d.sin2T4_bf==0.5 & d.mNu4Sq_bf>100);
Idxm4Max   = find( d.mNu4Sq_bf>1000);
fprintf('==========================================================================\n');
fprintf('%.0f  (%.0f significant)  best fits with sin2t4=0.5 and m4<100eV^2\n',numel(IdxSinMax1),sum(d.chi2_delta(IdxSinMax1)>5.99));
fprintf('%.0f   (%.0f significant)  best fits with sin2t4=0.5 and m4>100eV^2\n',numel(IdxSinMax2),sum(d.chi2_delta(IdxSinMax2)>5.99));
fprintf('%.0f (%.0f significant) best fits with m4>1000eV^2\n',numel(Idxm4Max),sum(d.chi2_delta(Idxm4Max)>5.99))
fprintf('==========================================================================\n');
% histogram(d.mNu4Sq_bf(Idx))