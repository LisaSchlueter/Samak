ST1=MultiRunAnalysis('DataType','Real','RunList','KITNUSCAN','exclDataStart',1,'fixPar','1 5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov','StackTolerance',0.2);
ST1.qUDistribution('saveplot','ON');

ST2=MultiRunAnalysis('DataType','Real','RunList','KITNUSCAN','exclDataStart',2,'fixPar','1 5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov','StackTolerance',0.1);
ST2.qUDistribution('saveplot','ON');


ST3=MultiRunAnalysis('DataType','Real','RunList','KITNUSCAN','exclDataStart',13,'fixPar','1 5 6','chi2','chi2Stat','Debug','ON','ELossFlag','Abdurashitov','StackTolerance',0.08);
ST3.qUDistribution('saveplot','ON');
