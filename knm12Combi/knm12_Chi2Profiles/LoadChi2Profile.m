function d = LoadChi2Profile(varargin)
p = inputParser;
p.addParameter('DataSet','Knm2',@(x)ismember(x,{'Knm1','Knm2'}));
p.addParameter('DataType','Real',@(x)ismember(x,{'Real','Twin'}));
p.addParameter('AnaStr','Uniform',@(x)ismember(x,{'Uniform','Ring_Full'}));
p.addParameter('chi2','chi2CMShape',@(x)ismember(x,{'chi2Stat','chi2CMShape'}));
p.addParameter('nFit',20,@(x)isfloat(x));
p.addParameter('mNuSqMin',-2.6,@(x)isfloat(x));
p.addParameter('mNuSqMax',1,@(x)isfloat(x));

p.parse(varargin{:});
DataSet  = p.Results.DataSet;
DataType = p.Results.DataType;
AnaStr   = p.Results.AnaStr;
chi2     = p.Results.chi2;
nFit     = p.Results.nFit;
mNuSqMin = p.Results.mNuSqMin;
mNuSqMax = p.Results.mNuSqMax;

filedir = [getenv('SamakPath'),'tritium-data/fit/',DataSet,'/Chi2Profile/',AnaStr,'/'];

switch DataSet
    case 'Knm1'
        NP = 1.064;
        if strcmp(chi2,'chi2Stat')
            chi2Str = chi2;
        else
            chi2Str = 'chi2CMShape_SysBudget24';
        end
    case 'Knm2'
        NP = 1.112;
        if strcmp(chi2,'chi2Stat')
            chi2Str = chi2;
        else
            if strcmp(AnaStr,'Uniform')
                chi2Str = 'chi2CMShape_SysBudget40';
            else
                chi2Str = 'chi2CMShape_SysBudget41';
            end
        end
end

filename = sprintf('%sChi2Profile_%s_UniformScan_mNu_%s_%sFPD_%s_NP%.3f_FitParE0BkgNorm_nFit%.0f_min-%.3g_max%.3g.mat',...
    filedir,DataType,DataSet,AnaStr,chi2Str,NP,nFit,abs(mNuSqMin),mNuSqMax);
try
    d = importdata(filename);
    fprintf('load from file: %s \n',filename)
catch
    fprintf(2,'file not found: %s \n',filename)
    d = 0;
    return
end

end