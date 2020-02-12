RunList = 'KNM1';
A = MultiRunAnalysis('RunList',RunList,'DataType','Twin','exclDataStart',2);
A.exclDataStart = 2;
A.Fit;
