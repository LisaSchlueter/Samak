%% Define figures positioning and basic default properties
close all;
%set(0,'DefaultFigureWindowStyle','normal');
set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultFigureColor','w');

RainbowFlag = true;

crsz = get(0,'ScreenSize');
figure_size_factor = 1;
std_figure_size = figure_size_factor*[640 410];
large_figure_size = figure_size_factor*[640 274];
if crsz(4)>900
    XY_scr = [1+1440/4 1+2*1024+1024/4];
else
    XY_scr = [1+1440/4 1+991/4];
end

set(0,'DefaultFigurePosition',[XY_scr std_figure_size]);

%% Initialization of the analysis
nobs = 40;
norm = 3000;

sys.t_acc=.1;
sys.t_li9=.03;
sys.t_prc=.1;
sys.s_acc=.01;
sys.s_li9 = .02;
sys.s_prc = .02;
sys.dm2 = 2.4e-3;
sys.s_dm2 = .05;

%% Fake data preparation
x=linspace(10/(2*nobs),10-10/(2*nobs),nobs)';
sigma_b2b = 0.01;

t_acc = sys.t_acc;
s_acc = sys.s_acc;
t_li9 = sys.t_li9;
s_li9 = sys.s_li9;
t_prc = sys.t_prc;
s_prc = sys.s_prc;
dm2 = sys.dm2;
s_dm2 = sys.s_dm2;
s2t13 = .1;

nu = (x>1).*(x-1).^2.*exp(-(x-1));
ampbds = (x>1).*(sin(1.27*dm2*(1+s_dm2*randn)*1.05e3./(x+0.782)).^2);
ampb = (x>1).*(sin(1.27*dm2*1.05e3./(x+0.782)).^2);
dm2_sys = -2*s2t13*(1.27*dm2*1.05e3./(x+0.782)).^2.*...
    sqrt(ampb.*(1-ampb));
nu = nu/sum(nu);
acc = (x>.5).*exp(-x);
acc = acc/sum(acc);
li9 = (x>1).*(x<10).*(x-1).*(10-x);
li9 = li9/sum(li9);
prc = (1+.2*randn(size(x))).*(x>.5).*ones(size(x));
prc = prc/sum(prc);

y=norm*(nu.*(1-s2t13*ampbds).*(1+1./(sqrt(norm*nu.*(1-s2t13*ampb))+eps).*...
    randn(size(x))+sigma_b2b*randn(size(x)))) + ...
    norm*(t_acc*acc.*(1+s_acc*randn(size(x))) + ...
    t_li9*li9.*(1+s_li9*randn(size(x))) + ...
    t_prc*prc.*(1+s_prc*randn(size(x))));

A=-norm*ampb.*nu;
S=norm*[nu acc li9 prc dm2_sys];
Winv = diag(1./[[.03 s_acc s_li9 s_prc] s_dm2].^2);
x0=[1 t_acc t_li9 t_prc dm2]';
Vinv=diag(1./(abs(y)+(norm*nu.*(1-s2t13*ampb)*sigma_b2b).^2+eps));

nsys = size(S,2);
npar = size(A,2);

%% Chi2 Analysis Tool with Systematics

s=cats(y,A,'Vinv',Vinv,'S',S,'Winv',Winv,'x0',x0);

%% Save interesting variables to base workspace

assignin('base','x',x);
assignin('base','y',y);
assignin('base','A',A);
assignin('base','s',s);
assignin('base','Vinv',Vinv);
assignin('base','Winv',Winv);
assignin('base','x0',x0);


%% Data and Model plot
figure(1);
clf;
set(gcf,'OuterPosition',[XY_scr large_figure_size]);
grid on;
subplot(1,5,[1 4]);
dataplot = [S(1:nobs,5)*s.xhat(5) S(1:nobs,4)*s.xhat(4) S(1:nobs,3)*s.xhat(3)...
    S(1:nobs,2)*s.xhat(2) S(1:nobs,1)*s.xhat(1)+A(1:nobs,1)*s.xhat(nsys+1)];
bar(x,dataplot,'stacked');
colormap(flipud(lines(5)));
hold on;
stairs([x; x(end)+(x(end)-x(end-1))]-10/(2*nobs),...
    [s.yhat(1:nobs); s.yhat(nobs)],'Color',.5*ones(3,1),'LineWidth',2);
stairs([x; x(end)+(x(end)-x(end-1))]-10/(2*nobs),...
    [s.yhat(1:nobs)-A(1:nobs,1)*s.xhat(nsys+1);...
    s.yhat(nobs)-A(nobs,1)*s.xhat(nsys+1)],'Color','k','LineWidth',2);
if(RainbowFlag)
    cmap = jet(nobs);
    for i=1:nobs,
        errorbar(x(i),y(i),1/sqrt(Vinv(i,i)),'ks','LineWidth',1,'MarkerFacecolor',cmap(i,:),...
        'MarkerSize',4);
    end
else
    errorbar(x,y,1./sqrt(diag(Vinv)),'ks','LineWidth',1,'MarkerFacecolor','y',...
        'MarkerSize',4);
end
hold off;
axis([0 10 min(min(min(dataplot,[],1))*1.2,0) max(y*1.2)]);
set(gca,'Xtick',0:10);
xlabel('E_{vis} in MeV');
ylabel('N_{evt}');
title('Data and Model');
legend('Dm2 Systematic','Proton Recoils','Li9','Accidentals','\nu signal',...
    'Total Oscillation Model','No oscillation Model','Data','Location',...
    'EastOutside');

subplot(1,5,5);
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
ini_p = [x0;s2t13]; % to check behaviour only with s2t13. In principle,
% no prior on its value...
pulls = zeros(nsys+npar,1);
for i=1:length(s.xhat)
    pulls(i) = (s.xhat(i)-ini_p(i))./sqrt(s.covx(i,i));
    bar(i,pulls(i),'FaceColor',cmap(i,:));
end
%errorbar(1:length(x0),x0,1./sqrt(diag(Winv)),'ko','LineWidth',2,'MarkerFac
%ecolor','y','MarkerSize',10);
axis([0 length(x0)+2 -1.2*max(abs(pulls)) 1.2*max(abs(pulls))]);
box on;
hold off;
title('Reduced Pulls');


%% Residuals and pulls plot
figure(2);
clf;
colormap(lines);
set(gcf,'OuterPosition',[XY_scr large_figure_size]);
subplot(1,3,1);
bar(s.r(1:nobs));
axis([0 length(s.r(1:nobs))+1 min(s.r(1:nobs))*1.2 max(s.r(1:nobs))*1.2]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.r(1:nobs)),10)));
xlabel('Bin number');
ylabel('R');
title('Residuals');

subplot(1,3,2);
qqplot(s.r(1:nobs));
box on;
grid on;
title('Residual QQ-plot');

subplot(1,3,3);
colormap(lines(5));
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
for i=1:nsys
    bar(i,s.r(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([0 length(s.r(nobs+(1:nsys)))+1 ...
    min(s.r(nobs+(1:nsys)))*1.2 max(s.r(nobs+(1:nsys)))*1.2]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.r(nobs+(1:nsys))),min(8,nsys))));
xlabel('Bin number');
ylabel('R');
title('Pulls');

%% Residuals vs Data
figure(3);
clf;
set(gcf,'OuterPosition',[XY_scr std_figure_size]);
hold on;
if(RainbowFlag)
    cmap = jet(nobs);
    resid = s.r(1:nobs);
    for i=1:nobs,
        plot(y(i),resid(i),'ko','LineWidth',1,'MarkerFacecolor',cmap(i,:),...
        'MarkerSize',8);
    end
else
    plot(y,s.r(1:nobs),'bo');
end
hold off;
box on;
% Think more about this plot...
%hold on; plot(S*x0,s.r(1:nobs),'rx'); hold off;
grid on;
xlabel('N_{evt}');
ylabel('R(y)');
title('Residuals vs Data');

%% Residuals vs Systematic variables
figure(4);clf;

ncol = ceil(sqrt(nsys));
nrow = ceil(nsys/ncol);

for i=1:nsys,
    %figure(3+i);
    
    subplot(nrow,ncol,i);

    hold on;

    sys_var = S(:,i)*s.xhat(i);
    
    if(RainbowFlag)
        cmap = jet(nobs);
        resid = s.r(1:nobs);
        for j=1:nobs,
            plot(sys_var(j),resid(j),'ko','LineWidth',1,'MarkerFacecolor',cmap(j,:),...
                'MarkerSize',6);
        end
    else
        plot(sys_var,resid,'bo');
    end
    hold off;
    box on;
    grid on;
    xlabel('N_{evt}');
    ylabel(['R(#' num2str(i) ')']);
    title(['Residuals vs Systematic var #' num2str(i)]);
end

%% Residuals vs Physical variables
figure(5);clf;

ncol = ceil(sqrt(npar));
nrow = ceil(npar/ncol);
resid = s.r(1:nobs);

for i=1:npar,
    %figure(3+nsys+i);
    
    subplot(nrow,ncol,i);
    
    hold on;

    phys_var = A(:,i)*s.xhat(nsys+i);
     
    if(RainbowFlag)
        cmap = jet(nobs);
        for j=1:nobs,
            plot(phys_var(j),resid(j),'ko','LineWidth',1,'MarkerFacecolor',cmap(j,:),...
                'MarkerSize',6);
        end
    else
        plot(phys_var,resid,'bo');
    end
    hold off;
    box on;
    grid on;
    xlabel('N_{evt}');
    ylabel(['R(#' num2str(i) ')']);
    title(['Residuals vs Physical var #' num2str(i)]);
end

%% Leverages plot
figure(4+nsys+npar);
clf;
bar(s.leverage(1:nobs));
colormap(lines(5));
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
for i=1:nsys
    bar(nobs+i,s.leverage(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([0 length(s.leverage)+1 0 max(s.leverage)*1.2]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.leverage),min(10,nobs+nsys))));
xlabel('E_{vis} in MeV');
ylabel('L');
title('Leverages');
line([1-1 length(s.leverage)+1],2*(nsys+npar)/(nobs+nsys)*[1 1],...
    'LineWidth',2,'LineStyle','--','Color',.5*ones(3,1));
line([1-1 length(s.leverage)+1],(nsys+npar)/(nobs+nsys)*[1 1],...
    'LineWidth',2,'LineStyle','-.','Color',.7*ones(3,1));

%% Delete-1 variance
figure(5+nsys+npar);
clf;
plot(s.s2_i,'kd','MarkerSize',10,'MarkerFaceColor','y','LineWidth',2);
axis([0 length(s.s2_i)+1 s.mse-6/sqrt(nobs-npar-1) s.mse+6/sqrt(nobs-npar-1)]);
line([1-1 length(x)+1],[s.mse s.mse],'LineWidth',2,'Color',.5*ones(3,1));
grid off;
set(gca,'Xtick',round(linspace(1,length(s.s2_i),min(8,nobs+nsys))));
xlabel('E_{vis} in MeV');
ylabel('\sigma_{i}');
title('Delete-1 variance');
if nobs+nsys-npar>5
    line([1-1 length(s.s2_i)+1],(s.mse+2/sqrt(nobs-npar-1))*[1 1],...
        'LineWidth',2,'LineStyle','--','Color',.5*ones(3,1));
    line([1-1 length(s.s2_i)+1],(s.mse-2/sqrt(nobs-npar-1))*[1 1],...
        'LineWidth',2,'LineStyle','--','Color',.5*ones(3,1));
end

%% Parameter correlation matrix
figure(6+nsys+npar);
clf;
nD=length(s.covx);
XD=1:nD;
YD=1:nD;
imagesc(XD,YD,sqrt(diag(1./diag(s.covx)))*s.covx*sqrt(diag(1./diag(s.covx))));
set(gca,'XTick',XD)
set(gca,'XTickLabel',{num2str(XD')})
set(gca,'YTick',YD)
set(gca,'YTickLabel',{num2str(YD')})
%set(gca,'YLabel',YD,YD-.5);
caxis([-1 1]);
cm=0:1/9:1;
my_cm_blue=[1-cm;1-cm;ones(size(cm))]';
my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
colorbar;
set(gcf,'Color',[1 1 1]);
% set(gca,'GridLineStyle','-'); set(gca,'MinorGridLineStyle','-');
% set(gca,'XGrid','on'); set(gca,'XMinorGrid','on'); set(gca,'YGrid','on');
% set(gca,'YMinorGrid','on');
set(gca,'XAxisLocation','top');
set(gca,'XTickLabel',{'nu' 'acc' 'li9' 'prc' 'dm2' 's2t13'})
set(gca,'YTickLabel',{'nu' 'acc' 'li9' 'prc' 'dm2' 's2t13'})
%grid on;
title('Parameter Correlation matrix');

%% Delete-1 parameters
figure(7+nsys+npar);
clf;
XD=1:size(s.x_i,2);
YD=1:size(s.x_i,1);
diffxi = s.x_i-repmat(s.xhat,[1 size(s.x_i,2)]);
imagesc(XD,YD,diffxi./repmat(max(abs(diffxi),[],2),[1 size(s.x_i,2)]));
set(gca,'Xtick',round(linspace(1,length(XD),min(8,nobs+nsys))));
set(gca,'YTick',round(linspace(1,length(YD),min(8,nsys+npar))));
caxis([-1 1]);
cm=0:1/9:1;
my_cm_blue=[1-cm;1-cm;ones(size(cm))]';
my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
colorbar;
set(gcf,'Color',[1 1 1]);
colorbar;
set(gca,'XAxisLocation','top');
grid on;
set(gca,'YTickLabel',{'nu' 'acc' 'li9' 'prc' 'dm2' 's2t13'})
title('Delete-1 parameters (centered and rescaled)');

%% Scale change in parameters
figure(8+nsys+npar);
clf;
XD=1:size(s.dfxs,2);
YD=1:size(s.dfxs,1);
imagesc(XD,YD,s.dfxs);
set(gca,'Xtick',round(linspace(1,length(XD),min(8,nobs+nsys))));
set(gca,'YTick',round(linspace(1,length(YD),min(8,nsys+npar))));
caxis([-1 1]);
cm=0:1/9:1;
my_cm_blue=[1-cm;1-cm;ones(size(cm))]';
my_cm_red=[ones(size(cm));1-cm;1-cm]';
colormap([flipud(my_cm_blue);[1 1 1];my_cm_red]);
colorbar;
set(gcf,'Color',[1 1 1]);
colorbar;
set(gca,'XAxisLocation','top');
grid on;
set(gca,'YTickLabel',{'nu' 'acc' 'li9' 'prc' 'dm2' 's2t13'})
title('Scale change in parameters');

%% Change in fitted values. DFFIT
figure(9+nsys+npar);
clf;
bar(s.dffit(1:nobs));
colormap(lines(5));
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
for i=1:nsys
    bar(nobs+i,s.dffit(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([0 length(s.dffit)+1 1.2*max(abs(s.dffit))*[-1 1]]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.dffit),min(10,nobs+nsys))));
xlabel('Bin number');
ylabel('dffit');
title('Change in fitted values. DFFIT');

%% Scaled change in fitted values. DFFITS (with studentized residuals)
figure(10+nsys+npar);
clf;
bar(s.dffits(1:nobs));
colormap(lines(5));
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
for i=1:nsys
    bar(nobs+i,s.dffits(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([0 length(s.dffits)+1 1.2*max(abs(s.dffits))*[-1 1]]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.dffits),min(10,nobs+nsys))));
xlabel('Bin number');
ylabel('dffits');
title('Scaled change in fitted values. DFFITS');

%% Change in covariance. COVRATIO
figure(11+nsys+npar);
clf;
bar(s.covratio(1:nobs),'b');
colormap(lines(5));
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
for i=1:nsys
    bar(nobs+i,s.covratio(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;

axis([0 length(s.covratio)+1 0 1.2*max(s.covratio)]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.covratio),min(10,nobs+nsys))));
xlabel('Bin number');
ylabel('covratio');
title('Change in covariance. COVRATIO');

%% Cook's Distance
figure(12+nsys+npar);
clf;
bar(s.cookd(1:nobs));
colormap(lines(5));
cmap = [.5*ones(1,3); circshift(lines(5),-1)];
hold on;
for i=1:nsys
    bar(nobs+i,s.cookd(nobs+i),'FaceColor',cmap(i,:));
end
hold off;
box on;
axis([0 length(s.cookd)+1 0 max(s.cookd)*1.2]);
grid off;
set(gca,'Xtick',round(linspace(1,length(s.cookd),min(10,nobs+nsys))));
xlabel('Bin number');
ylabel('Cookd');
title('Cook s Distance');
line([1 length(s.cookd)],[1 1],'LineWidth',2,'LineStyle','--',...
    'Color',.5*ones(3,1));

%% Residuals vs Leverages
figure(13+nsys+npar);
clf;
lev_edges = [0 min(max(s.leverage)*1.4,1)];
tmp_res_edges = max(abs(min(s.standres)),abs(max(s.standres)));
res_edges = [-tmp_res_edges tmp_res_edges]*1.4;
lev = linspace(lev_edges(1),lev_edges(2),500);
res = linspace(res_edges(1),res_edges(2),500);
[Lev,Res] = meshgrid(lev,res);
D = abs(Res).^2/(npar).*Lev./(1-Lev);
contourf(lev,res,D,[0 .5 1]);
line(2*(nsys+npar)/(nobs+nsys)*[1 1],res_edges,'LineWidth',2,...
    'LineStyle','-.','Color',.5*ones(3,1));
line((nsys+npar)/(nobs+nsys)*[1 1],res_edges,'LineWidth',2,...
    'LineStyle','--','Color',.7*ones(3,1));
line(lev_edges,[2 2],'LineWidth',2,'LineStyle','--','Color',.7*ones(3,1));
line(lev_edges,-[2 2],'LineWidth',2,'LineStyle','--','Color',.7*ones(3,1));
line(lev_edges,[3 3],'LineWidth',2,'LineStyle','-.','Color',.5*ones(3,1));
line(lev_edges,-[3 3],'LineWidth',2,'LineStyle','-.','Color',.5*ones(3,1));
axis([lev_edges res_edges]);
caxis([0 1]);
mycm=[[1 1 1];[8/9 8/9 1];[1 8/9 8/9]];
colormap(mycm);
%colorbar;
hold on;
if(RainbowFlag)
    cmap = jet(nobs);
    for i=1:nobs,
        plot(s.leverage(i),s.standres(i),'ko','LineWidth',1,'MarkerSize',10,...
            'MarkerFaceColor',cmap(i,:));
    end
    
    for i=1:nsys,
        cmap = [.5*ones(1,3); circshift(lines(5),-1)];
        plot(s.leverage(nobs+i),s.standres(nobs+i),'ks','LineWidth',1,'MarkerSize',10,...
            'MarkerFaceColor',cmap(i,:));
    end

else
    plot(s.leverage,s.standres,'kd','LineWidth',2,'MarkerSize',10,...
        'MarkerFaceColor','y');
end
grid on;
hold off;
xlabel('Leverages');
ylabel('Standardized Residuals');
title('Standardized Residuals vs Leverages');
hold on;
cookd_sel = s.cookd>=1;
cookd_sel_i = find(cookd_sel);
plot(s.leverage(cookd_sel),s.standres(cookd_sel),'ro','MarkerSize',20,...
    'LineWidth',2);
textx_offset = .04*(lev_edges(2)-lev_edges(1));
texty_offset = .02*(res_edges(2)-res_edges(1));
sys_txt = {'norm' 'acc' 'li9' 'prc' 'dm2'};
mysel = find(cookd_sel);
for i=1:length(mysel)
    if(mysel(i)<=nobs)
        mytxt = [num2str(cookd_sel_i(i)) ', E = ' sprintf('%1.1f',round(x(mysel(i))*10)/10) ' MeV'];
    else
        mytxt = sys_txt{mysel(i)-nobs};
    end
    text(s.leverage(cookd_sel_i(i))+textx_offset,...
        s.standres(cookd_sel_i(i))+texty_offset,...
        mytxt,'Color','r','FontWeight','bold');
end
hold off;

%% Correlation in Hat matrix
plotcorh(s)

%% T-statistics

TTable = dataset({s.tstat.xhat,'Coef'},{s.tstat.se,'StdErr'},...
    {s.tstat.t,'tStat'},{s.tstat.pval,'pVal'})
