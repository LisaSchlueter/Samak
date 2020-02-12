function stdrvslplot(s,RainbowFlag)
%% Residuals vs Leverages

if nargin<2
    RainbowFlag=true;
end

nobs = s.nobs;
nsys = s.nsys;
npar = s.npar;

clf;
lev_edges = [0 min(max(s.leverage)*1.4,1)];
tmp_res_edges = max(abs(min(s.standres)),abs(max(s.standres)));
res_edges = [-tmp_res_edges tmp_res_edges]*1.4;
lev = linspace(lev_edges(1),lev_edges(2),500);
res = linspace(res_edges(1),res_edges(2),500);
[Lev,Res] = meshgrid(lev,res);
D = abs(Res).^2/(npar+nsys).*Lev./(1-Lev);
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
% for i=1:length(mysel)
%     if(mysel(i)<=nobs)
%         mytxt = [num2str(cookd_sel_i(i)) ', E = ' sprintf('%1.1f',round(x(mysel(i))*10)/10) ' MeV'];
%     else
%         mytxt = sys_txt{mysel(i)-nobs};
%     end
%     text(s.leverage(cookd_sel_i(i))+textx_offset,...
%         s.standres(cookd_sel_i(i))+texty_offset,...
%         mytxt,'Color','r','FontWeight','bold');
% end
hold off;
PrettyFigureFormat;
end