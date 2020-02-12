% Compute Uncertainty on slope
function [CM CMFrac] = ComputeCM_Background(varargin)
p=inputParser;
p.addParameter('StudyObject','',isa('TBD'));
StudyObject = p.Results.SO;

TimeSec = 3*365*24*60*60;
qUfrac = 0.2;% time spent in background region
nqU = 3; % 3 Background points
qU  = [10 20 30]'+18575;
BKG_asimov = 0.3*TimeSec*qUfrac.*ones(nqU,1);
nSamples = 100;
par = zeros(2,nSamples);
err = zeros(2,nSamples);
chi2min = zeros(nSamples,1);

qULong = [-60:10:30]'+18575;
Bkg_Data = zeros(3,nSamples);
Bkg_Fit = zeros(numel(qULong),nSamples);

for i=1:nSamples
BKG = BKG_asimov + sqrt(BKG_asimov).*randn(nqU,1);
Data = [qU,BKG, sqrt(BKG)];
BKG_i = wmean(BKG,sqrt(BKG));
Slope_i = 0;
parInit = [BKG_i, Slope_i];
% Minuit Arguments
tmparg = sprintf(['set pri -10;'...
'migrad minos'],'');
% Minuit Input: Init Fit Parameter, Data, Covariance Matrix
Args = {parInit, Data, '-c', tmparg};
[par(:,i), err(:,i),chi2min(i), ~] = fminuit('chi2BKG',Args{:});
Bkg_Data(:,i) = Data(:,2)./(TimeSec*qUfrac); %for plot in cps
Bkg_Fit(:,i) =  ComputeBkgSlope(par(:,i),qULong)/(TimeSec*qUfrac);
end
%%

BkgIS = Bkg_Fit(:,chi2min<6);
plot(qULong-18575,1e3*BkgIS(:,:)); hold on;
plot(qU-18575,1e3.*Bkg_Data(:,1:100),'x'); 
PrettyFigureFormat;
xlabel('retarding potential - 18575 (eV)');
ylabel('Background (mcps)');
% legend(sprintf('BKG = %.2f \\pm %.2f mcps',1e3*par(1)/(TimeSec*qUfrac),1e3*err(1)/(TimeSec*qUfrac)),...
%     sprintf('Slope = %.2g \\pm %.2g Background Rate per eV',par(2)/(TimeSec*qUfrac),err(2)/(TimeSec*qUfrac)))
%%

CM = cov(BkgIS');
CMFrac = CM./0.3^2;
imagesc(CMFrac);
colorbar
%end


