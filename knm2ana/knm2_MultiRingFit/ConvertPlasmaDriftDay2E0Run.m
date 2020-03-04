function [E0,RectWidth] = ConvertPlasmaDriftDay2E0Run(varargin)
p = inputParser;
p.addParameter('DriftPerDay',[6*1e-03, 0, 6*1e-03]',@(x)all(isfloat(x)));
p.addParameter('E0OffseteV',[0,0.1,-0.1]',@(x)all(isfloat(x)));
p.addParameter('E0ref',18573.70,@(x)isfloat(x));
p.addParameter('SanityPlot','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

DriftPerDay = p.Results.DriftPerDay;  % in eV for 3 RW periods
E0OffseteV  = p.Results.E0OffseteV;   % step from RW-period to RW-period
E0ref       = p.Results.E0ref;        % average endpoint
SanityPlot  = p.Results.SanityPlot;

savedir = [getenv('SamakPath'),'knm2ana/knm2_MultiRingFit/results/'];
savename = sprintf('%sMCplasmaDriftE0_E0ref%.3feV_Drifts-%.0fmeV_%.0fmeV_%.0fmeV_Steps--%.0fmeV_%.0fmeV_%.0fmeV.mat',...
    savedir, E0ref,...
    DriftPerDay(1)*1e3,DriftPerDay(2)*1e3,DriftPerDay(3)*1e3,...
    E0OffseteV(1)*1e3,E0OffseteV(2)*1e3,E0OffseteV(2)*1e3);
if exist(savename,'file') 
    load(savename)
else
    MR = MultiRunAnalysis('RunList','KNM2_Prompt','DataType','Real');
    
    FirstRunRW2 = 56560;
    FirstRunRW3 = 57015;
    
    %% get duration
    StartTimeStamp   = MR.SingleRunData.StartTimeStamp;
    KNM2LiveTimeDays = days(StartTimeStamp-StartTimeStamp(1))+0.5*MR.SingleRunData.TimeSec./(60*60*24); % add half of measurement time to get middle time of run
    KNM2LiveTimeDaysRel = KNM2LiveTimeDays-mean(KNM2LiveTimeDays);
    
    %% convert Plasma Drift in Endpoint shift per run
    E0Shift  = -DriftPerDay.*KNM2LiveTimeDaysRel+E0OffseteV; % negative sign, because pos. plasma potential drift --> endpoint decreases
    
    E0 = zeros(MR.nRuns,1);
    E0(MR.RunList<FirstRunRW2) = E0ref + E0Shift(1,MR.RunList<FirstRunRW2);
    E0(MR.RunList>=FirstRunRW2 & MR.RunList<FirstRunRW3) = E0ref + E0Shift(2,MR.RunList>=FirstRunRW2 & MR.RunList<FirstRunRW3);
    E0(MR.RunList>=FirstRunRW3) = E0ref + E0Shift(3,MR.RunList>=FirstRunRW3);
    
    RectWidth = zeros(3,1);
    tmp1 = find(MR.RunList<FirstRunRW2); % end RW1
    tmp2 = find(MR.RunList>=FirstRunRW2 & MR.RunList<FirstRunRW3); % end RW2
    RectWidth(1) = E0(1)-E0(tmp1(end));
    RectWidth(2) = E0(tmp1(end)+1)-E0(tmp2(end));
    RectWidth(3) = E0(tmp2(end)+1)-E0(end);
    
    save(savename,'E0','KNM2LiveTimeDays','KNM2LiveTimeDaysRel','E0Shift','StartTimeStamp','FirstRunRW2','FirstRunRW3','RectWidth');
end
    if strcmp(SanityPlot,'ON')
        f2 = figure('Units','normalized','Position',[0.1,0.1,0.5,0.5]);
        plot(KNM2LiveTimeDays,E0-mean(E0),'s','MarkerSize',3,'MarkerFaceColor',rgb('DodgerBlue'),'Color',rgb('DodgerBlue'));
        PrettyFigureFormat('FontSize',24);
        xlabel('Live time (days)');
        ylabel(sprintf('E_0^{MC} - \\langleE_0^{MC}\\rangle  (eV)'));
        xlim([min(KNM2LiveTimeDays)-1 max(KNM2LiveTimeDays)+1]);
    end
end
