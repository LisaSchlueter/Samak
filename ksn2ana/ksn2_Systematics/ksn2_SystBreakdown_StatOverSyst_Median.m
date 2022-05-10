% ksn2 systematics breakdown
% median values for  ratio

DataType = 'Twin';
nGridSteps = 30;
range = 40;
InterpMode = 'spline';
CL = chi2cdf(1,1);

savedir = [getenv('SamakPath'),'ksn2ana/ksn2_Systematics/results/'];
savename = sprintf('%sksn2_SystBreakdown_StatOverSyst_%s_%.0feV_RasterScan%s_%sInterp_%.2gCL.mat',...
    savedir,DataType,range,'ON',InterpMode,CL);
if exist(savename,'file')
    d = importdata(savename);
    fprintf('load file %s \n',savename);
else
    fprintf('file not found %s \n Run ksn2_SystBreakdown_StatOverSyst.m',savename);
end

%%
LogIdx = d.mNu4Sq(1,:)<=600;
%% stat only
StatAll = median((d.sin2t4_Stat(1,:).^2./d.sin2t4_Tot(1,:).^2)');
StatLow = median((d.sin2t4_Stat(1,LogIdx).^2./d.sin2t4_Tot(1,LogIdx).^2)');
StatUp  = median((d.sin2t4_Stat(1,~LogIdx).^2./d.sin2t4_Tot(1,~LogIdx).^2)');
%% syst
SysAll = median(d.Ratio');
SysLow = median(d.Ratio(:,LogIdx)');
SysUp  = median(d.Ratio(:,~LogIdx)');

Summary = [StatAll,StatLow,StatUp;SysAll',SysLow',SysUp'];

%% print
fprintf('===================================================================================  \n');
fprintf('                                   Median all      ,    Median m4^2<600 eV^2   , Median m4^2>600 eV^2\n');
fprintf('Stat. only               :           %.2f         ,     %.2f                  , %.2f \n',Summary(1,:));
for i=1:d.nSys
    labeltmp = d.SysEffectLabel{i};
    if numel(labeltmp)<25
        labeltmp = [labeltmp,'                                    '];
    end
fprintf('%s:           %.2f         ,     %.2f                  , %.2f \n',labeltmp(1:25),Summary(i+1,:))
end
fprintf('===================================================================================  \n');
