function cor = cov2cor(cov)
SS = sqrt(diag(cov));
cor = cov./(SS*SS');
end
