%% -------------------------------------------------------------------------------------
%convert covariance matri into shape-only covariance matrix
% input:
% necessary:
% 1. CovMat, which should be converted into shape only
% 2. TBDIS from model, used for normalization:
% 3: background level
% optional: 
% exclDatastart to define energy range
% output: shape only covariance matrix and fractional shapy only covariance matrix
% Lisa, October 2019
%% -------------------------------------------------------------------------------------
function [CovMatShapeOnly, CovMatFracShapeOnly]= Convert2ShapeOnlyCovMat(CovMat,obj,varargin)
p=inputParser;
p.addParameter('exclDataStart','',@(x)isfloat(x));
p.parse(varargin{:});
exclDataStart = p.Results.exclDataStart;


TBDIS_NoBKG  = TBDIS-(BKG_RateSec.*obj.ModelObj.TimeSec.*obj.ModelObj.qUfrac);


end