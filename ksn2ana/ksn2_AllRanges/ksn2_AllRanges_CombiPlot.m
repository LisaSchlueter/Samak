
freePar = 'E0 Norm Bkg';
KSN1 = 'ON';
KSN2 = 'ON';
Mainz   = 'OFF';
Troitsk = 'OFF';
%% load knm2
savedir2 = [getenv('SamakPath'),'ksn2ana/ksn2_AllRanges/results/'];
savename90_knm2 = sprintf('%sksn2_AllRanges_chi2CMShape_%s_90eV.mat',savedir2,strrep(freePar,' ',''));
savename40_knm2 = sprintf('%sksn2_AllRanges_chi2CMShape_%s_40eV.mat',savedir2,strrep(freePar,' ',''));
fprintf('load file %s \n',savename90_knm2);
fprintf('load file %s \n',savename40_knm2);
dknm2_90 = importdata(savename90_knm2);
dknm2_40 = importdata(savename40_knm2);

%% load knm1
savedir1 = [getenv('SamakPath'),'ksn1ana/ksn1_ReAna/results/'];
savename90_knm1 = sprintf('%sksn1_AllRanges_chi2CMShape_%s_94eV.mat',savedir1,strrep(freePar,' ',''));
savename40_knm1 = sprintf('%sksn1_AllRanges_chi2CMShape_%s_40eV.mat',savedir1,strrep(freePar,' ',''));
fprintf('load file %s \n',savename90_knm1);
fprintf('load file %s \n',savename40_knm1);
dknm1_90 = importdata(savename90_knm1);
dknm1_40 = importdata(savename40_knm1);

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
    pMainz = plot(dMainz.SinSquare2Theta_X,dMainz.DmSquare41_Y,'-.','LineWidth',1.5,'Color',rgb('Black'));
    hold on;
end
% Troitsk
if strcmp(Troitsk,'ON')
    pTroitsk = plot(dTroitsk.SinSquare2Theta_X,dTroitsk.DmSquare41_Y,'--','LineWidth',1.5,'Color',rgb('DimGray'));
    hold on;
end

%%
if strcmp(KSN1,'ON')
    p1_40 = plot(dknm1_40.sin2T4_contour,dknm1_40.mNu4Sq_contour,':','LineWidth',3,'Color',rgb('Salmon'));
    hold on;
    p1_90_99 = plot(dknm1_90.sin2T4_contour99,dknm1_90.mNu4Sq_contour99,'-.','LineWidth',2.5,'Color',rgb('Silver'));
    
    if size(dknm1_90.sin2T4_contour,2)>1
        for i=1:size(dknm1_90.sin2T4_contour,2)
            p1_90 = plot(dknm1_90.sin2T4_contour(:,i),dknm1_90.mNu4Sq_contour(:,i),'-','LineWidth',2.5,'Color',rgb('Crimson'));
        end
    end
    
    
end

if strcmp(KSN2,'ON')
    p2_40 = plot(dknm2_40.sin2T4_contour,dknm2_40.mNu4Sq_contour,':','LineWidth',3,'Color',rgb('Turquoise'));
    hold on;
    p2_90 = plot(dknm2_90.sin2T4_contour,dknm2_90.mNu4Sq_contour,'-','LineWidth',3,'Color',rgb('MediumBlue'));
end

% if ~contains(freePar,'mNu')
%      savename = sprintf('%sksn12Combi_MaxRange_chi2CMShape_%s.mat',savedir2,strrep(freePar,' ',''));
%      dCombi90 = importdata(savename);
%      savedirCombi = sprintf('%sksn2ana/ksn2_RunCombination/results/',getenv('SamakPath'));
%     
%   %  savefileCombimNuFix = sprintf('%sksn21_Combination_ReAna_%s.mat',savedirCombi,'Real');
%   %  dCombi40 = importdata(savefileCombimNuFix);
%   %   p12_40 = plot(dCombi40.sin2T4_contour_12,dCombi40.mNu4Sq_contour_12,':','Color',rgb('ForestGreen'),'LineWidth',2);
%   %   p12_90 = plot(dCombi90.sin2T4_contour,dCombi90.mNu4Sq_contour,'-','Color',rgb('LimeGreen'),'LineWidth',2);
%   %   plot(dCombi90.BestFit.sin2T4,dCombi90.BestFit.mNu4Sq,'o','Color',rgb('LimeGreen'),'MarkerFaceColor',rgb('LimeGreen'),'MarkerSize',10);
% end

% Best fits
plot(dknm1_90.BestFit.sin2T4,dknm1_90.BestFit.mNu4Sq,'h','MarkerSize',15,'Color',p1_90.Color,'MarkerFaceColor',p1_90.Color);
%    plot(dknm1_40.BestFit.sin2T4,dknm1_40.BestFit.mNu4Sq,'o','MarkerSize',10,'Color',p1_40.Color,'MarkerFaceColor',p1_40.Color);
%  plot(dknm2_40.BestFit.sin2T4,dknm2_40.BestFit.mNu4Sq,'p','MarkerSize',15,'Color',p2_40.Color,'MarkerFaceColor',p2_40.Color);
plot(dknm2_90.BestFit.sin2T4,dknm2_90.BestFit.mNu4Sq,'p','MarkerSize',15,'Color',p2_90.Color,'MarkerFaceColor',p2_90.Color);


set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',22);

xlim([1e-03 0.5])
ylim([1 1e4]);
%%
pNone = plot(NaN,NaN,'.','Color','w');
pNone2 = plot(NaN,NaN,'.','Color','w');
pNone3 = plot(NaN,NaN,'.','Color','w');
pNone4 = plot(NaN,NaN,'.','Color','w');
%

%% 
leg = legend([pNone,pNone2,p1_40,p2_40,p1_90,p2_90,p1_90_99,pNone3],...
    'KNM1:','KNM2:',...
    '39 eV','40 eV',...
    '93 eV','90 eV',...
    '93 eV at 99% C.L.','',...
    'Location','southwest');
leg.NumColumns = 4;
PrettyLegendFormat(leg);
legend box off;
leg.Position(1)  = 0.1;


if contains(freePar,'mNu')
    tStr = sprintf('Analysis case II) {\\itm}_\\nu^2 free');
else
    tStr = sprintf('Analysis case I) {\\itm}_\\nu^2 = 0 eV^2');
end
t = text(1.2e-03,7.5,tStr,'FontSize',leg.FontSize);

%%
pltdir  = strrep(savedir2,'results','plots');
MakeDir(pltdir);
pltname = sprintf('%sContour_ExtendedRange_%s.pdf',pltdir,strrep(freePar,' ',''));
export_fig(pltname);
fprintf('save plot to %s \n',pltname);


