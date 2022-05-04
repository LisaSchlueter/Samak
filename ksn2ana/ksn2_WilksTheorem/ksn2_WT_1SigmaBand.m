% plot with 1 sigma sensitivity band 
% load randmoized mc contours

DataSet = 'KNM1';
if strcmp(DataSet,'KNM2')
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];
savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interplin_1000samples.mat',savedir);
else
    savedir = [getenv('SamakPath'),'ksn1ana/ksn1_WilksTheorem/results/'];
savefile = sprintf('%sksn1_WilksTheorem_40range_chi2CMShape_Mix.mat',savedir);
end
d = importdata(savefile);
% load ksn1 and 2 results
CombiDir = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));
fileT = sprintf('%sksn21_Combination_ReAna_Twin.mat',CombiDir);
fileD = sprintf('%sksn21_Combination_ReAna_Real.mat',CombiDir);
dCT = importdata(fileT);
dCD = importdata(fileD);

%% interpolate
mNu4Sq_contour = d.mNu4Sq_contour(~d.ClosedLog95);
sin2T4_contour = d.sin2T4_contour(~d.ClosedLog95);

mNu4Sq_min = max(cell2mat(cellfun(@(x) min(min(x)),mNu4Sq_contour,'UniformOutput',false))');
mNu4Sq_max = min(cell2mat(cellfun(@(x) max(max(x)),mNu4Sq_contour,'UniformOutput',false))');
mNu4Sq = logspace(log10(mNu4Sq_min),log10(mNu4Sq_max),1e3);
sin2T4 = zeros(sum(~d.ClosedLog95),1e3);

 
 InclIdx = true(sum(~d.ClosedLog95),1);
    for i=1:sum(~d.ClosedLog95)
        progressbar(i./sum(~d.ClosedLog95));
        x = mNu4Sq_contour{i};
        y = sin2T4_contour{i};
        if size(x,2)>1 && size(x,2)<500
            InclIdx(i) = 0;
        else
            sin2T4(i,:) = interp1(x,y,mNu4Sq,'spline','extrap');
        end
    end
    
    
    sin2T4 = sin2T4(InclIdx,:);
    
      GetFigure;
 sin2T4mean  = mean(sin2T4);
 sin2T4std   = std(sin2T4);
[l,a] = boundedline(sin2T4mean,mNu4Sq,sin2T4std,'orientation','horiz');
l.delete;
hold on;
pNone = plot(NaN,NaN,'w','LineStyle','none');
a.FaceColor = rgb('LightGray');
hold on;
if strcmp(DataSet,'KNM2')
    pD = plot(dCD.sin2T4_contour_2,dCD.mNu4Sq_contour_2,'-','LineWidth',3.5,'Color',rgb('DodgerBlue'));
    pT = plot(dCT.sin2T4_contour_2,dCT.mNu4Sq_contour_2,'-.','LineWidth',2.5,'Color',rgb('DimGray'));
else
    pD = plot(dCD.sin2T4_contour_1,dCD.mNu4Sq_contour_1,'-','LineWidth',3.5,'Color',rgb('SkyBlue'));
    pT = plot(dCT.sin2T4_contour_1,dCT.mNu4Sq_contour_1,'-.','LineWidth',2.5,'Color',rgb('DimGray'));
end
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([3e-03 0.5]);
ylim([mNu4Sq_min mNu4Sq_max])
PrettyFigureFormat('FontSize',22)
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
xlabel(sprintf('|{\\itU}_{e4}|^2'));
leg = legend([pNone,pT,a,pD],...
    sprintf('\nCase I) {\\itm}_\\nu^2 = 0 eV^2'),...
    'Sensitivity (Asimov)',...
    sprintf('1\\sigma sensitivity band'),...
    sprintf('%s data exclusion',DataSet),...
    'Location','southwest');PrettyLegendFormat(leg);


CombiPltDir = strrep(savedir,'results','plots');
MakeDir(CombiPltDir);
pltname = sprintf('%s%s_WT_1SigmaBand.png',CombiPltDir,DataSet);
print(pltname,'-dpng','-r350');
fprintf('save plot to %s \n',pltname);

