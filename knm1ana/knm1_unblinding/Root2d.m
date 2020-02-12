function F = Root2d(x,mNuSq,CumProb,DeltaChi2)
  F(1) =  (interp1(mNuSq,CumProb,x(2),'spline')-interp1(mNuSq,CumProb,x(1),'spline'))-0.9;
  F(2) =  interp1(mNuSq,DeltaChi2,x(1),'spline')-interp1(mNuSq,DeltaChi2,x(2),'spline');       
end