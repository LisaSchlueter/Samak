

DataType = 'Real';
freePar = 'E0 Bkg Norm';
range = 40;


savedir = sprintf('%sSterileAnalysis/GridSearchFiles/%s/%s/',...
    getenv('SamakPath'),'Knm1',DataType);

f1 = sprintf('%sKSN1_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid_Budget24%s.mat',...
    savedir,'KNM1',DataType,strrep(freePar,' ',''),range,'chi2CMShape',50,'_Negsin2T4');
f2 = sprintf('%sKSN1_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid_Budget24%s.mat',...
    savedir,'KNM1',DataType,strrep(freePar,' ',''),range,'chi2CMShape',50,'_NegmNu4Sq');
f3 = sprintf('%sKSN1_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid_Budget24%s.mat',...
    savedir,'KNM1',DataType,strrep(freePar,' ',''),range,'chi2CMShape',50,'_NegmNu4Sq_Negsin2T4');
f4 = sprintf('%sKSN1_GridSearch_%s_%s_%s_%.0feVrange_%s_%.0fnGrid_Budget24%s.mat',...
    savedir,'KNM1',DataType,strrep(freePar,' ',''),range,'chi2CMShape',50,'');
d1 = importdata(f1);
d2 = importdata(f2);
d3 = importdata(f3);
d4 = importdata(f4);


sin2T4 = [d1.sin2T4,d4.sin2T4;d3.sin2T4,d2.sin2T4];
mNu4Sq = [d1.mnu4Sq,d4.mnu4Sq;d3.mnu4Sq,d2.mnu4Sq];
chi2   = [d1.chi2',d4.chi2';d3.chi2',d2.chi2'];
chi2_ref = min(min(chi2));
DeltaChi2 = GetDeltaChi2(95,2);


%% find best fit per 
chi2tmp = d1.chi2';
[row, col]    = find(chi2tmp == min(chi2tmp(:)));
mNu4Sq_bf1 =   d1.mnu4Sq(row,col);%obj.mNu4Sq(col,row);
sin2T4_bf1 = d1.sin2T4(row,col);
chi2_bf1   =   min(min(chi2tmp));
            
chi2tmp = d2.chi2';
[row, col]    = find(chi2tmp == min(chi2tmp(:)));
mNu4Sq_bf2 =   d2.mnu4Sq(row,col);%obj.mNu4Sq(col,row);
sin2T4_bf2 = d2.sin2T4(row,col);
chi2_bf2  =   min(min(chi2tmp));

chi2tmp = d3.chi2';
[row, col]    = find(chi2tmp == min(chi2tmp(:)));
mNu4Sq_bf3 =   d3.mnu4Sq(row,col);%obj.mNu4Sq(col,row);
sin2T4_bf3 = d3.sin2T4(row,col);
chi2_bf3   =   min(min(chi2tmp));

chi2tmp = d4.chi2';
[row, col]    = find(chi2tmp == min(chi2tmp(:)));
mNu4Sq_bf4 =   d4.mnu4Sq(row,col);%obj.mNu4Sq(col,row);
sin2T4_bf4 = d4.sin2T4(row,col);
chi2_bf4   =   min(min(chi2tmp));
%% prepare plot

zlimMax = DeltaChi2;
 d1.chi2((d1.chi2-chi2_ref)>DeltaChi2) =  NaN;
 d2.chi2((d2.chi2-chi2_ref)>DeltaChi2) =  NaN;
d3.chi2((d3.chi2-chi2_ref)>DeltaChi2) =  NaN;
 d4.chi2((d4.chi2-chi2_ref)>DeltaChi2) =  NaN;

ymax = 2e3;

GetFigure;
s1 = subplot(2,2,1); % top left
surf(-d1.sin2T4,d1.mnu4Sq,d1.chi2'-chi2_ref,'EdgeColor','interp','FaceColor','interp');
set(gca,'XScale','log');
set(gca,'YScale','log');
view(180,-90); grid off;

PrettyFigureFormat('FontSize',20);
xlim([1e-03 0.5])
ylim([1 ymax])
xticks([0.001 0.01 0.1])
yticks([1 10 1e2 1e3]);
yticklabels({'','10^1','10^2','10^3'});
xticklabels({''});
ax1 = gca;
ylabel(sprintf('{\\itm}_4^2 (eV^2)'));
ctmp =colorbar;
caxis([0 zlimMax])
ctmp.delete;

s2 = subplot(2,2,2); % top right
surf(d4.sin2T4,d4.mnu4Sq,d4.chi2'-chi2_ref,'EdgeColor','interp','FaceColor','interp');
view(2)
set(gca,'XScale','log');
set(gca,'YScale','log');
xlim([0 0.5])
PrettyFigureFormat('FontSize',20);
grid off;
xlim([1e-03 0.6])
ylim([1 ymax])
yticks([1 10 1e2 1e3]);
xticks([0.001 0.01 0.1])
yticklabels('');
xticklabels('');
ax2 = gca;
ctmp =colorbar;
caxis([0 zlimMax])
ctmp.delete;
hold on;
if strcmp(DataType,'Real')
    pbf = plot3(sin2T4_bf4,mNu4Sq_bf4,10,'x','MarkerSize',9,'Color',rgb('White'),'LineWidth',2);
    leg = legend(pbf,'\color{white}Best fit', 'EdgeColor','none','Location','southwest','Color',rgb('White'),...
        'FontSize',get(gca,'FontSize'));
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;0.3]));
end
s3 = subplot(2,2,4); % bottom right
surf(d2.sin2T4,-d2.mnu4Sq,d2.chi2'-chi2_ref,'EdgeColor','interp','FaceColor','interp');
view(2);
set(gca,'XScale','log');
 set(gca,'YScale','log');
view(0,-90);
xlim([1e-03 0.6]);
ylim([1 ymax])
PrettyFigureFormat('FontSize',20);
grid off;
yticks([1 10 1e2 1e3]);
xticks([0.001 0.01 0.1])
yticklabels({''});
ax3 = gca;
c =colorbar;
c.Label.String = sprintf('\\Delta\\chi^2');
c.Label.FontSize = get(gca,'FontSize')+2;
c.Limits=[0 zlimMax];

s4 = subplot(2,2,3); % bottom left
surf(-d3.sin2T4,-d3.mnu4Sq,d3.chi2'-chi2_ref,'EdgeColor','interp','FaceColor','interp');
view(180,90);
set(gca,'XScale','log');
 set(gca,'YScale','log');
PrettyFigureFormat('FontSize',20)
grid off;
xlim([1e-03 0.5])
ylim([1 ymax])
xticks([0.001 0.01 0.1])
yticks([1 10 1e2 1e3]);
yticklabels({'-10^0','-10^1','-10^2','-10^3'});
xticklabels({'','-10^{-2}','-10^{-1}'});
ax4 = gca;
xlabel(sprintf('|{\\itU}_{e4}|^2'));
ctmp =colorbar;
caxis([0 zlimMax])
ctmp.delete;
%% move closer together
ax1.Position(2) = 0.61;
ax4.Position(2) = 0.245;
ax3.Position(2) = 0.245;
ax2.Position(2) = 0.61;
ax3.Position(1) = 0.48;
ax2.Position(1) = 0.48;
ax3.Position(3) = ax1.Position(3);

ax1.YLabel.Position(2) = 1;
ax1.YLabel.Position(1) = 1.5;
ax1.YLabel.FontSize = 21;
ax4.XLabel.Position(1) = 0.8e-03;
ax4.XLabel.FontSize = 21;

c.Position(1) = 0.83;
c.Position(2) = 0.3;
c.Position(4) = 0.6;
%%
plotdir = sprintf('%sSterileAnalysis/plots/%s/%s/', getenv('SamakPath'),'Knm1',DataType);
MakeDir(plotdir);
plotname = sprintf('%sNegParSpace_%s_%.0feV_%s',plotdir,DataType,range,strrep(freePar,' ',''));
print(gcf,[plotname,'.png'],'-dpng','-r300');
%export_fig(gcf,[plotname,'.pdf']);
print(gcf,[plotname,'.png'],'-dpng','-r300');
fprintf('save plot to %s \n',plotname);
