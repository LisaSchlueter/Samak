%
nGrid = 50;
nSamples = 1e3;
Mode = 'NH';
savedir = [getenv('Samak3.0'),'ksn1ana/ksn1_NuBetaBeta/results/'];

%% load files
range = 40;
DataType = 'Real';
freePar = 'mNu E0 Norm Bkg';
pullFlag = 15;
switch pullFlag
    case 99
        pullStr = '';
    case 15
        pullStr = sprintf('_mNuSqPull1eV2');
end

file1 = sprintf('%sContour_Real_%s_40eV.mat',savedir,'E0NormBkg');
file2 = sprintf('%sContour_Real_%s_40eV.mat',savedir,'mNuE0NormBkg');
file3 = sprintf('%sContour_Real_%s_40eV_mNuSqPull1eV2.mat',savedir,'mNuE0NormBkg');
file4 = sprintf('%sContour_FinalSensitivity.mat',savedir);
d1 = importdata(file1);
d2 = importdata(file2);
d3 = importdata(file3);
d4 = importdata(file4);
%%
savename = sprintf('%sksn1_NuBetaBeta_ToyMC_%s_Grid%.0f_Samples%.0f.mat',savedir,Mode,nGrid,nSamples);

if exist(savename,'file') %&&1==2
    load(savename)
else

% grid

m4Sq    = repmat(logspace(-1,4,nGrid),nGrid,1);
sint4Sq = repmat(logspace(-3,log(0.74),nGrid),nGrid,1)';

phase    = 2.*rand(1,1,nSamples)-1;  % draw random phase for each grid point between +-1

if strcmp(Mode,'NH')
    % normal hierarchy: possible values in non-degenerate regime: 0-0.005 eV
    mbb3nu    = 0.005.*rand(1,1,nSamples);
elseif strcmp(Mode,'IH')
    % inverted hierarchy: possible values in non-degenerate regime: 0.01-0.05 eV
    mbb3nu    = 0.04.*rand(1,1,nSamples)+0.01;
end

mbbexp = abs((1-sint4Sq).*mbb3nu+exp(1i.*pi.*phase).*sqrt(m4Sq).*sint4Sq);

sint2Sq_contour =repmat(logspace(-3,log(0.74),nGrid*50),nSamples,1)';%zeros(nGrid*10,nSamples);%cell(nSamples,1);
m4Sq_contour   = zeros(nGrid*50,nSamples);
Index = ones(nSamples,1);
for i=1:nSamples
    if all(all(mbbexp(:,:,i)<0.2)) || all(all(mbbexp(:,:,i)>0.2))
        Index(i) = 0;
        continue
    end
    [a,~] = contour(sint4Sq,m4Sq,squeeze(mbbexp(:,:,i)),[0.2,0.2]);
    hold on;
    try
        m4Sq_contour(:,i) = interp1(a(1,2:end),a(2,2:end),sint2Sq_contour(:,i),'lin','extrap');
    catch
        Index(i) = 0;
    end
end
close all
%% 
Index = logical(Index);
m4Sq_contour= m4Sq_contour(:,Index);
sint2Sq_contour = sint2Sq_contour(:,Index);
%x = cell2mat(cellfun(@(x) x(2:end),sint2Sq_contour,'UniformOutput',false)');
%y = cell2mat(cellfun(@(x) x(2:end),m4Sq_contour,'UniformOutput',false)');%cell2mat(m4Sq_contour');
% %%
% % % calculate density for each sint2Sq_contour grid point
% % for i=1:nGrid*10
% %     sint2Sq_bin = sint2Sq_contour(i,1);
% %     
% %     sort(m4Sq_contour(i,:));
% %     
% % end
% x = reshape(sint2Sq_contour,numel(sint2Sq_contour),1);
% y = reshape(m4Sq_contour,numel(m4Sq_contour),1);
%  h1 = hexscatter(x,y,'o','MarkerEdgeColor','none','MarkerFaceColor',rgb('DodgerBlue'),'MarkerFaceAlpha',0.01,...
%      'SizeData',15);
% % close all
% % for i=1:nSamples/2
% % tmpx = sint2Sq_contour{i};
% % tmpy = m4Sq_contour{i};
% % a = area(tmpx(2:end),tmpy(2:end),1e3,'FaceAlpha',0.2,'EdgeColor','none','FaceColor',rgb('SkyBlue'));
%  set(gca,'YScale','log');
%  set(gca,'XScale','log');
% % hold on;
% % end
save(savename,'sint2Sq_contour','m4Sq_contour')
end
%%
p = plot(sint2Sq_contour,m4Sq_contour,'Color',[rgb('DodgerBlue'),0.05]);
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
PrettyFigureFormat('FontSize',20);
hold on;
pcontour    = plot(d1.sinT4Sq,d1.mNu4Sq,'-','Color',rgb('Orange'),'LineWidth',2.5);
pSensitivity = plot(d4.sinT4Sq,d4.mNu4Sq,':','Color',rgb('Silver'),'LineWidth',2.5);
pNone = plot(ones(10,1),NaN.*ones(10,1),'-','Color',rgb('DodgerBlue'));
xlim([1e-03 0.5]);
ylim([1 2e3]);
%%
if strcmp(Mode,'NH')
    legStr  = sprintf('{\\itm}_{\\beta\\beta}^{3\\nu} \\in [0 eV, 0.005 eV]');
else
    legStr  = sprintf('{\\itm}_{\\beta\\beta}^{3\\nu} \\in [0.01 eV, 0.05 eV]');
end
leg = legend([pNone,pcontour,pSensitivity], ...
    sprintf('{\\itm}_{\\beta\\beta}^{exp} = 0.2 eV'),...
    sprintf('KATRIN {\\itm}_\\nu^2 fixed (95%% C.L.)'),...
    sprintf('Projected KATRIN Final sensitivity {\\itm}_\\nu^2 fixed (95%% C.L.)'),...
    'Location','southwest');
leg.EdgeColor = rgb('Silver');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.7])); % make legend semi-transparent

leg.FontSize = 14;

title(sprintf('Randomly drawn \\delta \\in [-1,1] and %s',legStr),...
    'FontWeight','normal','FontSize',15);
%%
print(sprintf('./plots/ksn1_NuBetaBeta_ToyMC_%s.png',Mode),'-dpng','-r300');