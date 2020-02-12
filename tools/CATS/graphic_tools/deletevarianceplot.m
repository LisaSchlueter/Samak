function deletevarianceplot(s)
%% Delete-1 variance

plot(s.s2_i,'kd','MarkerSize',10,'MarkerFaceColor','y','LineWidth',2);
axis([0 length(s.s2_i)+1 s.mse-6/sqrt(s.nobs-s.npar-1) s.mse+6/sqrt(s.nobs-s.npar-1)]);
line([1-1 s.nobs+s.nsys+1],[s.mse s.mse],'LineWidth',2,'Color',.5*ones(3,1));
grid off;
set(gca,'Xtick',round(linspace(1,length(s.s2_i),min(8,s.nobs+s.nsys))));
xlabel('E_{vis} in MeV');
ylabel('\sigma_{i}');
title('Delete-1 variance');
if s.nobs+s.nsys-s.npar>5
    line([1-1 length(s.s2_i)+1],(s.mse+2/sqrt(s.nobs-s.npar-1))*[1 1],...
        'LineWidth',2,'LineStyle','--','Color',.5*ones(3,1));
    line([1-1 length(s.s2_i)+1],(s.mse-2/sqrt(s.nobs-s.npar-1))*[1 1],...
        'LineWidth',2,'LineStyle','--','Color',.5*ones(3,1));
end

end