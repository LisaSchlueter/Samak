Twin =MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',1,'fixPar','5 6 7 8 9 10 11');
Twin.Fit;
%%
Twin2 =MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',1,'fixPar','5 6 7 8 9 10 11','StackTolerance',0.1);
Twin2.Fit;

Twin3 =MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',1,'fixPar','5 6 7 8 9 10 11','StackTolerance',0.05);
Twin3.Fit;