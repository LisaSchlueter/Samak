% Test of Wilk's theorem (coverage)
% chi2 distribution of best fits

Hypothesis = 'H0';
MergeNew = 'ON';
RmDuplicates = 'ON';
InterpMode = 'Mix';
switch Hypothesis
    case 'H0'
        NrandMC = 419;
        Twin_sin2T4 = 0;
        Twin_mNu4Sq = 0;
        randMC_new  = 1:1250;
        MergeStr = sprintf('_MergeNew%.0f',numel(randMC_new));
    case 'H1'
        randMC = [1:1500];
        % excl = [1:139,577:757];
        % randMC = randMC(~ismember(randMC,excl));
        Twin_sin2T4 = 0.0240;
        Twin_mNu4Sq = 92.7;
        chi2 = 'chi2CMShape';
        MergeNew = 'OFF'; % nothing new
        NrandMC = numel(randMC);
end

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_WilksTheorem/results/'];

if Twin_sin2T4==0 && Twin_mNu4Sq==0
    savefile = sprintf('%sksn2_WilksTheorem_NullHypothesis_Interp%s_%.0fsamples%s_RmDouble%s.mat',...
        savedir,InterpMode,NrandMC,MergeStr,RmDuplicates);
else
    savefile = sprintf('%sksn2_WilksTheorem_mNu4Sq-%.1feV2_sin2T4-%.3g_Interp%s_%.0fsamples.mat',...
        savedir,Twin_mNu4Sq,Twin_sin2T4,InterpMode,NrandMC);
end

  savefileContour = strrep(savefile,'.mat','_Contour.mat');
  
if exist(savefileContour,'file')
    load(savefileContour)
else
    
    if exist(savefile,'file')
        load(savefile);
        fprintf('load file from %s \n',savefile);
    else
        fprintf('file does not exist: %s \n',savefile);
        return
    end
    
    
    
    %% interpolate contour at same mNu4Sq
    mNu4Sq_min = max(cell2mat(cellfun(@(x) min(min(x)),mNu4Sq_contour(~ClosedLog95),'UniformOutput',false))');
    mNu4Sq_max = min(cell2mat(cellfun(@(x) max(max(x)),mNu4Sq_contour(~ClosedLog95),'UniformOutput',false))');
    mNu4Sq = logspace(log10(mNu4Sq_min),log10(mNu4Sq_max),1e3);
    sin2T4 = zeros(sum(~ClosedLog95),1e3);
    
    Idx = find(~ClosedLog95);
    
    InclIdx = true(numel(Idx),1);
    for i=1:sum(~ClosedLog95)
        progressbar(i./numel(Idx));
        x = mNu4Sq_contour{Idx(i)};
        y = sin2T4_contour{Idx(i)};
        if size(x,2)>1 && size(x,2)<500
            InclIdx(i) = 0;
        else
            sin2T4(i,:) = interp1(x,y,mNu4Sq,'spline','extrap');
        end
    end
    
    sin2T4 = sin2T4(InclIdx,:);
    save(savefileContour,'sin2T4','mNu4Sq','sin2T4_contour_Asimov','mNu4Sq_contour_Asimov','mNu4Sq_min','mNu4Sq_max')
end
%%
GetFigure;
[l,a] = boundedline(mean(sin2T4),mNu4Sq,std(sin2T4),'orientation','horiz');
l.delete;
a.FaceColor = rgb('LightGray');
hold on;
pA = plot(sin2T4_contour_Asimov,mNu4Sq_contour_Asimov,'-','LineWidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([3e-03 0.5]);
ylim([mNu4Sq_min mNu4Sq_max])
PrettyFigureFormat('FontSize',22)
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
xlabel(sprintf('|{\\itU}_{e4}|^2'));
leg = legend([pA,a],'Asimov sensitivity',sprintf('Randomized MC: 1\\sigma band'),'Location','southwest');
PrettyLegendFormat(leg);
t = title(sprintf('Null hypothesis , {\\itm}_\\nu^2 = 0 eV^2'));
t.FontWeight = 'normal';t.FontSize = get(gca,'FontSize');
pltdir = strrep(savedir,'results','plots');
pltname = strrep(strrep(savefileContour,'results','plots'),'.mat','.png');
print(pltname,'-dpng','-r350');
%% load ksn-2
savedir = [getenv('SamakPath'),'ksn2ana/ksn2_DataTwin/results/'];
MakeDir(savedir);
savefile2 = sprintf('%sksn2_DataTwinContour_%s_%s.mat',savedir,'chi2CMShape','E0NormBkg');
if exist(savefile2,'file')
    d2 = importdata(savefile2);
    fprintf('load from file %s \n',savefile2);
    pD = plot(d2.sin2T4_contourD2,d2.mNu4Sq_contourD2,':','LineWidth',2,'Color',rgb('Red'));
    leg = legend([pA,a,pD],...
        'Asimov sensitivity',sprintf('Randomized MC: 1\\sigma band'),...
        'Data exclusion', 'Location','southwest');
    PrettyLegendFormat(leg);
    pltdir = strrep(savedir,'results','plots');
    pltname = strrep(strrep(savefileContour,'results','plots'),'.mat','_Data.png');
    print(pltname,'-dpng','-r350');
end

