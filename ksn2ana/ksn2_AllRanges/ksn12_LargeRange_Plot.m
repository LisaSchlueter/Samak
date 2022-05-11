
freePar = 'E0 Norm Bkg';
Mainz   = 'ON';
Troitsk = 'ON';
%% load knm2
savedir2 = [getenv('SamakPath'),'ksn2ana/ksn2_AllRanges/results/'];
savename90_knm2 = sprintf('%sksn2_AllRanges_chi2CMShape_%s_90eV.mat',savedir2,strrep(freePar,' ',''));
fprintf('load file %s \n',savename90_knm2);
d2 = importdata(savename90_knm2);

%% load knm1
savedir1 = [getenv('SamakPath'),'ksn1ana/ksn1_ReAna/results/'];
savename90_knm1 = sprintf('%sksn1_AllRanges_chi2CMShape_%s_94eV.mat',savedir1,strrep(freePar,' ',''));
fprintf('load file %s \n',savename90_knm1);
d1 = importdata(savename90_knm1);

%% load Maiinz / Troitsk
savedirOther = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];
% load Mainz
filenameMainz = sprintf('%scoord_Mainz_95CL.mat',savedirOther);
dMainz = importdata(filenameMainz);

% load Troitsk
filenameTroitsk = sprintf('%scoord_Troitsk_95CL.mat',savedirOther);
dTroitsk = importdata(filenameTroitsk);

%%
f1 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.45]);
if strcmp(Mainz,'ON')
     sinTsq = 0.5*(1-sqrt(1-dMainz.SinSquare2Theta_X));
    pMainz = plot(sinTsq,dMainz.DmSquare41_Y,'-.','LineWidth',1.5,'Color',rgb('Black'));
    hold on;
end
% Troitsk
if strcmp(Troitsk,'ON')
    pTroitsk = plot(dTroitsk.SinSquareTheta_X,dTroitsk.DmSquare41_Y,'--','LineWidth',1.5,'Color',rgb('DimGray'));
    hold on;
end

%%

%p1_99 = plot(d1.sin2T4_contour99,d1.mNu4Sq_contour99,'-.','LineWidth',2.5,'Color',rgb('Silver'));

if size(d1.sin2T4_contour,2)>1
    for i=1:size(d1.sin2T4_contour,2)
        p1 = plot(d1.sin2T4_contour(:,i),d1.mNu4Sq_contour(:,i),'-','LineWidth',2.5,'Color',rgb('Crimson'));
        hold on;
    end
end

p2 = plot(d2.sin2T4_contour,d2.mNu4Sq_contour,'-','LineWidth',2.5,'Color',rgb('MediumBlue'));



savename = sprintf('%sksn12Combi_MaxRange_chi2CMShape_%s.mat',savedir2,strrep(freePar,' ',''));
dCombi90 = importdata(savename);
savedirCombi = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));

p12 = plot(dCombi90.sin2T4_contour,dCombi90.mNu4Sq_contour,'-','Color',rgb('LimeGreen'),'LineWidth',3);


% Best fits
plot(dCombi90.BestFit.sin2T4,dCombi90.BestFit.mNu4Sq,'o','Color',rgb('LimeGreen'),'MarkerFaceColor',rgb('LimeGreen'),'MarkerSize',10);
plot(d1.BestFit.sin2T4,d1.BestFit.mNu4Sq,'h','MarkerSize',15,'Color',p1.Color,'MarkerFaceColor',p1.Color);
plot(d2.BestFit.sin2T4,d2.BestFit.mNu4Sq,'p','MarkerSize',15,'Color',p2.Color,'MarkerFaceColor',p2.Color);


set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);

xlim([1e-03 0.5])
ylim([1 1e4]);

leg = legend([p1,p2,p12(1),pMainz,pTroitsk],...
    'KNM1','KNM2','KNM1+2',...
'Mainz','Troitsk',...
    'Location','southwest');
leg.NumColumns = 2;
PrettyLegendFormat(leg);
legend box off;


%%
tStr = sprintf('Analysis case I) {\\itm}_\\nu^2 = 0 eV^2');
t = text(1.2e-03,15,tStr,'FontSize',leg.FontSize);

%
pltdir  = strrep(savedir2,'results','plots');
MakeDir(pltdir);
pltname = sprintf('%sContour_ksn12Combi_ExtendedRange_%s.pdf',pltdir,strrep(freePar,' ',''));
export_fig(pltname);
fprintf('save plot to %s \n',pltname);


