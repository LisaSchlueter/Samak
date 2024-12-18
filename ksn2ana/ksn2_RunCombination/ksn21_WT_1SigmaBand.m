% 1 sigma cnotour band from randomized mc study
% ksn1+2 combination
InterpMode = 'lin';
Twin_sin2T4 = 0;
Twin_mNu4Sq = 0;
RecomputeFlag = 'OFF';
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn21_WilksTheorem_NullHypothesis_Interp%s.mat',...
        savedir,InterpMode);
else
    savefile = sprintf('%sksn21_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode);
end

if exist(savefile,'file') && strcmp(RecomputeFlag,'OFF')
    fprintf('savefile already created \n');
   d = importdata(savefile);
  % load(savefile)
else
    return
end

%% load ksn12 result
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
    %%
    GetFigure;
 %%   ksn1+2 combined
 sin2T4mean  = mean(sin2T4);
 sin2T4std   = std(sin2T4);
[l,a] = boundedline(sin2T4mean,mNu4Sq,sin2T4std,'orientation','horiz');
l.delete;
hold on;
pNone = plot(NaN,NaN,'w','LineStyle','none');
a.FaceColor = rgb('LightGray');
hold on;
%pA = plot(mean(sin2T4),mNu4Sq,'-','LineWidth',2);
pD = plot(dCD.sin2T4_contour_12,dCD.mNu4Sq_contour_12,'-','LineWidth',3.5,'Color',rgb('Navy'));
pT = plot(dCT.sin2T4_contour_12,dCT.mNu4Sq_contour_12,'-.','LineWidth',2.5,'Color',rgb('DimGray'));
 pDbf = plot(dCD.sin2T4_bf_12,dCD.mNu4Sq_bf_12, 'p','MarkerIndices',[1,1],'MarkerFaceColor',pD.Color,...
       'LineWidth',1,'Color',pD.Color,'MarkerSize',17);
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
    'KNM1+2 data exclusion',...
    'Location','southwest');PrettyLegendFormat(leg);

CombiPltDir = strrep(CombiDir,'results','plots');
MakeDir(CombiPltDir);
pltname = sprintf('%sksn21_1SigmaBand.png',CombiPltDir);
print(pltname,'-dpng','-r350');
fprintf('save plot to %s \n',pltname);

