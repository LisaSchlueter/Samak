% RingWise Fit Using Ring Analysis Class

MyRingList=1:12;
ReDoFit='ON';
RunList = 'KNM1';

savepath = [getenv('SamakPath'),'tritium-data/fit/'];
if ~exist(savepath,'dir')
    system(['mkdir ',savepath]);
end

savefile = sprintf('BkgPerPix_%s_%.0fnPix',RunList,numel());

if ~exist('R','var')
        M = MultiRunAnalysis('RunList',RunList,'RingList',MyRingList);
        R = RingAnalysis('RunAnaObj',M,'RingList',MyRingList);
end

% Fit Rings
R.FitRings;

%%
BkgFit = (R.FitResult.par(:,3)+cell2mat(arrayfun(@(x) x.ModelObj.BKG_RateSec_i,R.MultiObj,'UniformOutput',0))'); % background from fit per ring
nPixFit = cell2mat(arrayfun(@(x) numel(x.PixList),R.MultiObj,'UniformOutput',0))'; %number of pixels per ring in fit

BkgPix = zeros(148,1); % relative Bkg per pixel

[~,RingPixList] = Ring2Pixel(MyRingList,1:148);
for i=1:numel(MyRingList)
    BkgPix(M.RingPixList{i}) = repmat(BkgFit(i),numel(M.RingPixList{i}),1)./nPixFit(i); % background per pixel
end

BkgRelPix = BkgPix./sum(BkgPix);



