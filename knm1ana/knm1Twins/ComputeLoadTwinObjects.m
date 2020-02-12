function [Twin, TwinqU, TwinqUfrac, TwinqUqUfrac, TwinIs, TwinCD, Twinall]...
    = ComputeLoadTwinObjects(varargin)
% function to calculate twins objects with different same flags
% with be saved in results

p=inputParser;
p.addParameter('RunList','KNM1',@(x)ischar(x));
p.parse(varargin{:});
RunList       = p.Results.RunList;

%% calculate/load twin objects
savedir = [getenv('SamakPath'),'knm1ana/knm1Twins/results/'];
if ~exist(savedir,'dir')
    system(['mkdir ',savedir]);
end

savename = ['Twins_',RunList];
if exist([savedir,savename,'.mat'],'file')
    fprintf(2,'Loading Twin Object \n')
    load([savedir,savename,'.mat']);
else
    Twin         = MultiRunAnalysis('RunList',RunList,'DataType','Twin','exclDataStart',17);
    save([savedir,savename,'.mat'],'Twin');
end

if exist([savedir,savename,'_qU.mat'],'file')
    fprintf(2,'Loading Twin Object - qU same \n')
    load([savedir,savename,'_qU.mat']);
else
    TwinqU       = MultiRunAnalysis('RunList',RunList,'DataType','Twin','Twin_SameqUFlag','ON','exclDataStart',17);
    save([savedir,savename,'_qU.mat'],'TwinqU');
end

if exist([savedir,savename,'_qUfrac.mat'],'file')
    fprintf(2,'Loading Twin Object - qUfrac same\n')
    load([savedir,savename,'_qUfrac.mat']);
else
    TwinqUfrac   = MultiRunAnalysis('RunList',RunList,'DataType','Twin','Twin_SameqUfracFlag','ON','exclDataStart',17);
    save([savedir,savename,'_qUfrac.mat'],'TwinqUfrac');
end

if exist([savedir,savename,'_qUqUfrac.mat'],'file')
    fprintf(2,'Loading Twin Object - qU+qUfrac same\n')
    load([savedir,savename,'_qUqUfrac.mat']);
else
    TwinqUqUfrac = MultiRunAnalysis('RunList',RunList,'DataType','Twin','Twin_SameqUFlag','ON','Twin_SameqUfracFlag','ON','exclDataStart',17);
    save([savedir,savename,'_qUqUfrac.mat'],'TwinqUqUfrac');
end

if exist([savedir,savename,'_CD.mat'],'file')
    fprintf(2,'Loading Twin Object rhod same \n')
    load([savedir,savename,'_CD.mat']);
else
    TwinCD      = MultiRunAnalysis('RunList',RunList,'DataType','Twin','Twin_SameCDFlag','ON','exclDataStart',17);
    save([savedir,savename,'_CD.mat'],'TwinCD');
end

if exist([savedir,savename,'_Is.mat'],'file')
    fprintf(2,'Loading Twin Object - isotropologues same \n')
    load([savedir,savename,'_Is.mat']);
else
    TwinIs        = MultiRunAnalysis('RunList',RunList,'DataType','Twin','Twin_SameIsotopFlag','ON','exclDataStart',17);
    save([savedir,savename,'_Is.mat'],'TwinIs');
end


if exist([savedir,savename,'_all.mat'],'file')
    fprintf(2,'Loading Twin Object - all same \n')
    load([savedir,savename,'_all.mat']);
else
    Twinall      = MultiRunAnalysis('RunList',RunList,'DataType','Twin','exclDataStart',17,...
        'Twin_SameCDFlag','ON','Twin_SameIsotopFlag','ON','Twin_SameqUFlag','ON','Twin_SameqUfracFlag','ON');
    save([savedir,savename,'_all.mat'],'Twinall');
end



end
