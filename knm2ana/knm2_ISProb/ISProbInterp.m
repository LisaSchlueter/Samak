

[RhoDSigma,Theta,ISProb] = Compute_InitISProbMeshGrid;

[X,Y] = meshgrid(RhoDSigma,Theta);
ISProb0 = squeeze(ISProb(1,:,:))';

surf(X,Y,ISProb0);
xlabel('rhodsigma');
ylabel('theta');
zlabel('Is prob 0');


IsProb0Inter = interp2(X,Y,ISProb0,mean(RhoDSigma),mean(Theta),'spline');
