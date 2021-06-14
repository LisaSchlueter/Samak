% contour with randomized twins
%% settings
range = 40;%
chi2Str = 'chi2CMShape';
InterpMode = 'Mix';
savedir = [getenv('SamakParh'),'ksn1ana/ksn1_WilksTheorem/results/'];
savename = sprintf('%sksn1_WilksTheorem_%.0frange_%s_%s.mat',savedir,range,chi2Str,InterpMode);
load(savename);



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
%% load ksn1 (Re-Analysis)
k1file = [getenv('SamakPath'),'ksn2ana/ksn2_DataTwin/results/ksn1_DataTwinContour_chi2CMShape_E0NormBkg.mat'];
d1 = importdata(k1file);

%%
GetFigure;
[l,a] = boundedline(mean(sin2T4),mNu4Sq,std(sin2T4),'orientation','horiz');
l.delete;
a.FaceColor = rgb('LightGray');
hold on;
%pA = plot(mean(sin2T4),mNu4Sq,'-','LineWidth',2);
pT = plot(d1.sin2T4_contourT1,d1.mNu4Sq_contourT1,'-','LineWidth',2);
pD = plot(d1.sin2T4_contourD1,d1.mNu4Sq_contourD1,'-.','LineWidth',2,'Color',rgb('FireBrick'));
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([3e-03 0.5]);
ylim([mNu4Sq_min mNu4Sq_max])
PrettyFigureFormat('FontSize',22)
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
xlabel(sprintf('|{\\itU}_{e4}|^2'));
leg = legend([pT,a,pD],'Asimov sensitivity',...
    sprintf('Randomized MC: 1\\sigma band'),...
    'Data exclusion',...
    'Location','southwest');
PrettyLegendFormat(leg);