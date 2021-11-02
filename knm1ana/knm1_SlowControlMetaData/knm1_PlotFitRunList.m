% knm1 original settings
% load single run fit results via object (uniform, fixed nu-mass)


savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    return
end

CommonArg = {'HideGaps','OFF','ShowRWPeriods','OFF','saveplot','OFF','TimeStyle','date'};
%R.PlotFitRunList(CommonArg{:},'Parameter','RhoD','YLim',1e17.*[1.076,1.15]);
R.PlotFitRunList('Parameter','T2');
%R.PlotFitRunList('HideGaps','OFF','ShowRWPeriods','OFF','YLim',[18572.6,18575.3],'Parameter','E0','saveplot','pdf');
%R.PlotFitRunList('HideGaps','OFF','ShowRWPeriods','OFF','Parameter','B','saveplot','pdf');
%R.PlotFitRunList('HideGaps','OFF','ShowRWPeriods','OFF','YLim',1+[-1,+1].*0.12,'Parameter','N','saveplot','pdf');