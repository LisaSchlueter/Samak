try
    M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11');
    S = RunSensitivity('RunAnaObj',M);
    S.ComputeFCR('nSamples',1000,'mNuSq_t',[0.1,0.2,0.3:0.2:0.9,5:0.5:6]);
catch
end
% 
% try [Twin, TwinqU, TwinqUfrac, TwinqUqUfrac, TwinIs, TwinCD, Twinall]...
%     = ComputeLoadTwinObjects(varargin);
% catch
% end
% 
% try
%    FitTwins('exclDataStart',14);
% catch
% end
% 
% try
%    FitTwins('exclDataStart',17);
% catch
% end
% 
% try
%    FitTwins('exclDataStart',2);
% catch
% end
% 
% try
%     M = MultiRunAnalysis('RunList','KNM1','DataType','Twin','exclDataStart',14,'fixPar','5 6 7 8 9 10 11');
%     S = RunSensitivity('RunAnaObj',M);
%     S.ComputeFCR('nSamples',1000,'mNuSq',[0.1:0.2:3],'FC','ON');
% catch
% end
