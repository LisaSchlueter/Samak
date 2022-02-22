function stdrlevplot(s,RainbowFlag,NameFlag)
% Residual versus leverage plot of the fit performed by CATS to identify
% indluential points.
%
% stdrvslplot(s,RainbowFlag,NameFlag) plots the studentized residuals versus
% the leverages of the fit performed by CATS where s is the output result of CATS
% ( s = cats(...) ). The Cook's distance is computed for each point and and
% Cook's distance level cruve are plotted (D = 0.5 and D = 1). By default, RainbowFlag is 
% set to true and colors the points from blue to red to identify begining to end of the
% spectrum in this plot as well as systematic bins. By default, NameFlag is
% set to true and indicates names of bins where the Cook's distance is
% above 1.
%
% Copyright 2005-2012 Guillaume MENTION, CEA Saclay
% $Revision: 1.1$  $Date: 2012/02/08$


nobs = s.nobs;
nsys = s.nsys;
npar = s.npar;

if nargin<2
    RainbowFlag = true;
end

clf;
lev_edges = [0 min(max(s.leverage)*1.4,1)];
tmp_res_edges = max(abs(min(s.standres)),abs(max(s.standres)));
res_edges = [-tmp_res_edges tmp_res_edges]*1.4;
lev = linspace(lev_edges(1),lev_edges(2),500);
res = linspace(res_edges(1),res_edges(2),500);
[Lev,Res] = meshgrid(lev,res);
D = abs(Res).^2/(npar+nsys).*Lev./(1-Lev);
contourf(lev,res,D,[0 .5 1]);
% line(2*(nsys+npar)/(nobs+nsys)*[1 1],res_edges,'LineWidth',2,...
%     'LineStyle','-.','Color',.5*ones(3,1));
% line((nsys+npar)/(nobs+nsys)*[1 1],res_edges,'LineWidth',2,...
%     'LineStyle','--','Color',.7*ones(3,1));
line(lev_edges,[0 0],'LineWidth',2,'LineStyle','-','Color',.7*ones(3,1));
% line(lev_edges,[2 2],'LineWidth',2,'LineStyle','--','Color',.7*ones(3,1));
% line(lev_edges,-[2 2],'LineWidth',2,'LineStyle','--','Color',.7*ones(3,1));
% line(lev_edges,[3 3],'LineWidth',2,'LineStyle','-.','Color',.5*ones(3,1));
% line(lev_edges,-[3 3],'LineWidth',2,'LineStyle','-.','Color',.5*ones(3,1));
axis([lev_edges res_edges]);
caxis([0 1]);
mycm=[[1 1 1];[8/9 8/9 1];[1 8/9 8/9]];
colormap(mycm);
%colorbar;
hold on;
if(RainbowFlag)
    cmap = jet(nobs);
    for i=1:nobs
        p(i) = plot(s.leverage(i),s.standres(i),'ko','LineWidth',1,'MarkerSize',14,...
            'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','none');
    end
    
    cmap = lines(nsys);
    for i=1:nsys
        plot(s.leverage(nobs+i),s.standres(nobs+i),'kp','LineWidth',1,'MarkerSize',14,...
            'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor','none');
    end

else
    plot(s.leverage,s.standres,'kd','LineWidth',2,'MarkerSize',14,...
        'MarkerFaceColor','y');
end 
grid on;
hold off;
xlabel('Leverages');
ylabel(sprintf('Standardized Residuals (\\sigma)'));
%title('Standardized Residuals versus Leverages (Cook''s Distance)','FontSize',18);
hold on;
cookd_sel = s.cookd>=1;
cookd_sel_i = find(cookd_sel);
 plot(s.leverage(cookd_sel),s.standres(cookd_sel),'ro','MarkerSize',20,...
    'LineWidth',2);

textx_offset = .04*(lev_edges(2)-lev_edges(1));
texty_offset = .02*(res_edges(2)-res_edges(1));
mysel = find(cookd_sel);
if ~isempty(s.names)
    for i=1:length(mysel)
        text(s.leverage(cookd_sel_i(i))+textx_offset,...
            s.standres(cookd_sel_i(i))+texty_offset,...
            s.names{cookd_sel_i(i)},'Color','r','FontWeight','bold');
    end
else
    for i=1:length(mysel)
        text(s.leverage(cookd_sel_i(i))+textx_offset,...
            s.standres(cookd_sel_i(i))+texty_offset,...
            sprintf('%d',cookd_sel_i(i)),'Color','r','FontWeight','bold');
    end    
end
hold off;
PrettyFigureFormat('FontSize',24);
%grid off
counter = 1;
leg_entry=cell(numel(p),1);
%leg_entry{1} = 'qU_{min}';
for i=1:numel(p)%2:numel(p)-1
    leg_entry{counter}=sprintf('%.0f eV',s.names(i));%sprintf('%0.f',counter);
    counter = counter +1;  
end
%leg_entry{end} = 'qU_{max}';
leg =legend(p,leg_entry{:});
leg.Location = 'northeastoutside';
leg.Title.FontWeight = 'normal';
leg.Title.String = sprintf('{\\itqU} - 18574 eV');
leg.NumColumns = 2;
legend boxoff;
end