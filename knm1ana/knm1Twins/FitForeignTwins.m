DataType = 'KafitTwin';%'FitriumTwin';%'KafitTwin';%'FitriumTwin';%'KafitTwin';%
exclDataStart = 14;
RunArg = {'RunList','KNM1','DataType',DataType,'fixPar','5 6 7 8 9 10 11',...
    'exclDataStart',exclDataStart};

if ~exist('M','var')
    M = MultiRunAnalysis(RunArg{:});
end