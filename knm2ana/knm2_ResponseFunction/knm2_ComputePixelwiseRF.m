%% compute pixelwiese response function for KNM2
% compute pixelwise spectrum
    % settings
RunList               = 'KNM2_Prompt';
exclDataStart         = 11; % 40eV range = 27 subruns
chi2                  = 'chi2Stat';
fixPar = 'E0 Bkg Norm';
pullFlag = 4;
% Init Model Object and covariance matrix object
savedir = [getenv('SamakPath'),'knm2ana/knm2_ResponseFunction/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_ResponseFunction_Pixelwise2.mat')];
Arg = {'RunList',RunList,...
    'chi2','chi2Stat','DataType','Real',...
    'exclDataStart',exclDataStart,...
    'fixPar',fixPar,...% free Parameter !
    'RadiativeFlag','ON',...
    'NonPoissonScaleFactor',1,...
    'minuitOpt','min ; minos',...
    'FSDFlag','BlindingKNM2',...
    'ELossFlag','KatrinT2',...
    'SysBudget',22,...
    'AnaFlag','StackPixel',...
    'chi2',chi2,...
    'pullFlag',pullFlag};
Rall = MultiRunAnalysis(Arg{:});
TimeSec_i = Rall.ModelObj.TimeSec;
PixList = Rall.PixList;

RF    = zeros(Rall.ModelObj.nTe,Rall.ModelObj.nqU,numel(PixList));
Model = cell(numel(PixList),1);
Te    = cell(numel(PixList),1);
qU    = zeros(Rall.ModelObj.nqU,numel(PixList));
TBDIS = zeros(Rall.ModelObj.nqU,numel(PixList));

progressbar('computing piselwise response function')
for i=1:numel(PixList)
    progressbar(i/numel(PixList));
    Model{i} = MultiRunAnalysis(Arg{:},'PixList',PixList(i));
    Model{i}.ModelObj.ComputeTBDDS;
    Model{i}.ModelObj.ComputeTBDIS;
    
    TBDIS(:,i) = Model{i}.ModelObj.TBDIS;
    Te{i}      = Model{i}.ModelObj.Te;
    qU(:,i)    = Model{i}.ModelObj.qU;
    RF(:,:,i)  = Model{i}.ModelObj.RF;
end
save(savename,'RF','Arg','PixList','Model','Te','qU','TBDIS','-mat');

%% stack pixelwise spectrum
TBDIS_Uniform = sum(TBDIS,2);
qU_Uniform    = mean(qU,2); % mean same as weighted mean here, because all pixels have the same time
save(savename,'TBDIS_Uniform','qU_Uniform','-append');
