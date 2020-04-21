function chi2 = Chi2GaussKNM12MultiRing(p,X)
% Gaussian Xi2
% X = [ "column vector of bin centers"  | "column vector of data" | "column vector of uncertainties" | "... extra info" ]
% p = [ p1; p2; ... ] column vector of parameter to fit

% Th. Lasserre - CEA Saclay - December 2012

%knm1 - ring loop
j=0;
for i=[0 9 18 27]
    j=j+1;
    knm1x(:,j) = X(i+1:i+7,1); knm2x(:,j) = X(8+i:9+i,1);
    knm1y(:,j) = X(i+1:i+7,2); knm2y(:,j) = X(8+i:9+i,2);
    knm1z(:,j) = X(i+1:i+7,3); knm2z(:,j) = X(8+i:9+i,3);
end

k1r1 = Model(p(1:2),knm1x(:,1));
k2r1 = Model(p(3:4),knm2x(:,1));

k1r2 = Model(p(5:6),knm1x(:,2));
k2r2 = Model(p(7:8),knm2x(:,2));

k1r3 = Model(p(9:10),knm1x(:,3));
k2r3 = Model(p(11:12),knm2x(:,3));

k1r4 = Model(p(13:14),knm1x(:,4));
k2r4 = Model(p(15:16),knm2x(:,4));


slopeErrorIntraKNM12 = 1e-2;
slopeErrorIntraPSR   = 1e-2;

chi2 = sum(( knm1y(:,1) - k1r1 ).^2 ./ knm1z(:,1).^2) + ...
       sum(( knm2y(:,1) - k2r1 ).^2 ./ knm2z(:,1).^2) + ...
       sum(( knm1y(:,2) - k1r2 ).^2 ./ knm1z(:,2).^2) + ...
       sum(( knm2y(:,2) - k2r2 ).^2 ./ knm2z(:,2).^2) + ...
       sum(( knm1y(:,3) - k1r3 ).^2 ./ knm1z(:,3).^2) + ...
       sum(( knm2y(:,3) - k2r3 ).^2 ./ knm2z(:,3).^2) + ...
       sum(( knm1y(:,4) - k1r4 ).^2 ./ knm1z(:,4).^2) + ...
       sum(( knm2y(:,4) - k2r4 ).^2 ./ knm2z(:,4).^2) + ...
       ((p(4)/p(3)-p(2)/p(1))/slopeErrorIntraKNM12).^2 + ...
       ((p(8)/p(7)-p(6)/p(5))/slopeErrorIntraKNM12).^2 + ...
       ((p(12)/p(11)-p(10)/p(9))/slopeErrorIntraKNM12).^2 + ...
       ((p(16)/p(15)-p(14)/p(13))/slopeErrorIntraKNM12).^2 + ...
       ((p(6)/p(5)-p(2)/p(1))/slopeErrorIntraPSR).^2 + ...
       ((p(10)/p(9)-p(6)/p(5))/slopeErrorIntraPSR).^2 + ...
       ((p(14)/p(13)-p(10)/p(9))/slopeErrorIntraPSR).^2;
end
