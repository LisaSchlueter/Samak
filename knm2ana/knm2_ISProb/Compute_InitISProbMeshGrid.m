function [RhoDSigma,Theta,ISProb] = Compute_InitISProbMeshGrid(varargin)
p=inputParser;
p.addParameter('DataSet','Knm2Theta+',@(x)ismember(x,{'Knm2','Knm1','Knm2Theta+','Knm1Theta+'}));
p.addParameter('NIS',7,@(x)isfloat(x));
p.parse(varargin{:});
DataSet = p.Results.DataSet;
NIS     = p.Results.NIS;

savedir = [getenv('SamakPath'),'inputs/WGTSMACE/WGTS_ISProb/'];
savename = [savedir,sprintf('InitISProbMeshGrid_%s.mat',DataSet)];

if NIS~=7
    savename = strrep(savename,'.mat',sprintf('_%.0fNIS.mat',NIS));
end

if exist(savename,'file')
    load(savename,'RhoDSigma','Theta','ISProb')
else
    
    switch DataSet
        case {'Knm1','Knm1Theta+'}
            A = ref_FakeRun_KNM1_CD22_23days('NIS',NIS);
        case {'Knm2','Knm2Theta+'}
            A = ref_FakeRun_KNM2_CD84_2hours('NIS',NIS);
    end
    %%
    maxEis = 1000;
    IsProbBinStep = 2;
    Eiscs = 18575+(-maxEis:IsProbBinStep:maxEis);
    if contains(DataSet,'Theta+')
        %         Bmax_Min = A.MACE_Bmax_T.*0.90;
        %         Bmax_Max = A.MACE_Bmax_T.*2;
        %         BmaxSamples = linspace(Bmax_Min,Bmax_Max,100);
        BmaxFun = @(theta) A.WGTS_B_T./sin(theta).^2;
        ThetaSamples = linspace(1e-05,1,100);
        BmaxSamples = BmaxFun(ThetaSamples);
    else
        Bmax_Min = A.MACE_Bmax_T.*0.90;
        Bmax_Max = A.MACE_Bmax_T.*1.10;
        BmaxSamples = linspace(Bmax_Min,Bmax_Max,20);
    end
    ISProb = zeros(A.NIS+1,numel(Eiscs),numel(BmaxSamples));
    
    for i=1:numel(BmaxSamples)
        progressbar(i/numel(BmaxSamples));
        file_pis =  sprintf('%sIS_%.5g-molPercm2_Edep-Xsection-max%.0feV_Xstep%.1feV_%.0f-NIS_%.3g-Bmax_%.3g-Bs.mat',...
            savedir,A.WGTS_CD_MolPerCm2,maxEis,IsProbBinStep,A.NIS+1,BmaxSamples(i),A.WGTS_B_T);
        if ~exist(file_pis,'file')
            A.MACE_Bmax_T = BmaxSamples(i);
            A.ComputeISProb('Energy',reshape(Eiscs,[1,1,numel(Eiscs)]),...
                            'Method','Exact');
        end
        
        d = importdata(file_pis);
        ISProb(:,:,i) = d.Pis_m(:,:);
    end
    
    RhoDSigma = A.WGTS_CD_MolPerCm2.*A.ISXsection(Eiscs);
    ThetaFun = @(bmax,bs)  asin(sqrt(bs./bmax));
    Theta    = ThetaFun(BmaxSamples,A.WGTS_B_T);
    save(savename,'RhoDSigma','Theta','ISProb');
end
end