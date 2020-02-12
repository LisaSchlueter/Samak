function [mNu_Fit, parScan, errScan, chi2minScan, dof ,mNu90,mNumin] = NuMassScan_SensitivityNominal_Nplus1Systematics_Loop1(varargin)
p = inputParser;
p.addParameter('SysBudget','03',@(x)ischar(x));
p.addParameter('MACE_Ba_T',7*1e-04,@(x)isfloat(x));
p.addParameter('WGTS_B_T',0.7*3.6,@(x)isfloat(x));
p.addParameter('range',30,@(x)isfloat(x));
p.addParameter('TimeSec',3*365*24*60*60,@(x)isfloat(x));
p.addParameter('RecomputeFlag','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('TD','',@(x)ischar(x)); %if empty, fill accordning to range and Ba ("optimal")
p.addParameter('ScanPrcsn',0.02,@(x)isfloat(x));
p.addParameter('mNuStop',0.6,@(x)isfloat(x));
p.addParameter('Scan','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('BKG_RateSec','',@(x)isfloat(x)); %if empty, filled according to Ba
p.addParameter('Q_i',18575,@(x)isfloat(x));
p.addParameter('PlotFit','OFF',@(x)ischar(x));
p.addParameter('SaveResults','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

SysBudget     = p.Results.SysBudget;
MACE_Ba_T     = p.Results.MACE_Ba_T;
WGTS_B_T      = p.Results.WGTS_B_T;
range         = p.Results.range;
TimeSec       = p.Results.TimeSec;
RecomputeFlag = p.Results.RecomputeFlag;
TD            = p.Results.TD;
ScanPrcsn     = p.Results.ScanPrcsn;
mNuStop       = p.Results.mNuStop;
Scan          = p.Results.Scan;
BKG_RateSec   = p.Results.BKG_RateSec;
Q_i           = p.Results.Q_i;
PlotFit       = p.Results.PlotFit;
SaveResults   = p.Results.SaveResults;
                           
mySysEffects = {'TC',...%TC only
    struct('TC','ON','FSD','ON'),...%TC + FSD
    'all'}; %TC + RF + FSD = everything
nFitMax = 20; 
    
if isempty(TD)
TD = sprintf('Sensitivity_%.0feV_Ba%.0fG',range,MACE_Ba_T*1e4);
end
TDlabel = strrep(TD,sprintf('_Ba%.0fG',MACE_Ba_T*1e4),''); %
TDlabel = strrep(TDlabel,'Sensitivity_','');
save_name = sprintf('./results/%s_SensitivityNominal_ResultsNuMassScan_MTD-%s_Ba%.0fG_Bs%.2fT_SystematicsNplus1.mat',SysBudget,TDlabel,MACE_Ba_T*1e4,WGTS_B_T);
   
    %Init: gather results
    chi2minScan  = zeros(numel(mySysEffects)+1,nFitMax); % chi2min distribution
    mNu90        = zeros(numel(mySysEffects)+1,1);                  % sensitivity on mnu^2 (90% C.L.)
    mNumin       = zeros(numel(mySysEffects)+1,1);                  % mass with minimal chi2
    parScan      = zeros(numel(mySysEffects)+1,6,nFitMax);
    errScan      = zeros(numel(mySysEffects)+1,6,nFitMax);
    mnuSq_i_Fit  = zeros(numel(mySysEffects)+1,nFitMax);
    %Load File if possible
    if exist(save_name,'file')==2 && strcmp(RecomputeFlag,'OFF')
        sprintf('File already exist. Do you want to recompute? \n')
        return
    else
        %% Do Scan
        if isempty(BKG_RateSec)
            BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T);
        end
        Scan_Arg = {'TimeSec',TimeSec,...
            'MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Scan',Scan,'BKG_RateSec',BKG_RateSec,...
            'range',range,'ScanPrcsn',ScanPrcsn,'mNuStop',mNuStop,'plotFit',PlotFit,'TD',TD,'Q_i',Q_i};
        %stat
        [mnuSq_i_Fit(1,:), parScan(1,:,:), errScan(1,:,:), chi2minScan(1,:,:),~,mNu90(1),mNumin(1)] = ...
            NuMassScan_SensitivityNominal(ScanArg{:},'chi2','chi2Stat');
        
        % Systematics
        for i=1:numel(mySysEffects)
            [mnuSq_i_Fit(i+1,:), parScan(i+1,:,:), errScan(i+1,:,:), chi2minScan(i+1,:,:), dof ,mNu90(i+1),mNumin(i+1)] = ...
                NuMassScan_SensitivityNominal(Scan_Arg{:},'chi2','chi2CM',...
                'SysEffect',mySysEffects{i},'SysBudget',SysBudget,'PlotCM','OFF');
        end
        % save       
        if strcmp(SaveResults,'ON')
            save(save_name,'parScan','errScan','chi2minScan','ScanPrcsn','mNu90', 'mNumin','TD','mySysEffects','dof','WGTS_B_T','MACE_Ba_T_all','BKG_RateSec','TimeSec');
        end
    end
end