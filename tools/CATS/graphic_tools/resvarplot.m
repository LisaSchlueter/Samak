function resvarplot(s)
%% Residuals vs Systematic and Physical variables
figure(4);clf;

ncol = ceil(sqrt(nsys));
nrow = ceil(nsys/ncol);
nsys = s.nsys;
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
end
