function CorrMat = ComputeCorrMatqUScan(exclnPoints,nqUmax)
      %% Calculate correlated mean and p-value
      % first: covariance matrix, because it's the same data set! Accordig
      % to Cowan book (page 112)
      
      % nqUmax: maximal number of subruns considered
      
      % qU Scan starts from excluding few subruns to exclude many subruns
       nCommon = zeros(exclnPoints);
       %nSingle = 12:-1:1;
       
       % number of considered subruns 
       nSingle = nqUmax:-1:(nqUmax-exclnPoints+1);
       
       for i=1:exclnPoints
           for j=1:exclnPoints
               if i>=j
                   nCommon(i,j) = nqUmax-i+1;
               else
                   nCommon(i,j) = nqUmax-j+1;
               end
               
           end
       end  
       
       CorrMat = nCommon./sqrt((nSingle'.*nSingle));
end      