% fit background level as a function of time
function knm2_FitBKGSlopeTime(obj)

switch obj.DataSet
    case 'Knm2'
        obj.exclDataStart = 11;
        obj.fixPar = 'E0 Bkg Norm';
        obj.InitFitPar;
end

if isempty(obj.SingleRun_FitResults.chi2Stat)
    obj.FitRunList('Recompute','ON')
end
end