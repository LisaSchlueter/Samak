function [qUslope  StartTimeStampMean] = GetMosDriftKNM2(varargin)
% drift of main spectrometer voltage as a function of time
% measurement with Krypton line position in monitor spectrometer (MoS)
% analysis done by Thibaut & Vikas 2020
% output of this script: main spectrometer qU drift in eV/day for KNM2

p=inputParser;
p.addParameter('SanityPlot','ON',@(x) ismember(x,{'ON','OFF'}));
p.parse(varargin{:});
SanityPlot = p.Results.SanityPlot;

savedir = [getenv('SamakPath'),'inputs/MonitorSpec/'];
savename = sprintf('%sKnm2_Mos_qUDrift.mat',savedir);

if exist(savename,'file')
    load(savename);
else
    d = importdata([savedir,'linepos_mos.txt']);
    
    RunNr     = d(:,1);
    Time      = (d(:,2)-d(1,2))./(60*60*24); % days
    KrPos     = d(:,3);
    KrPosErr  = d(:,4);
    
    
    % cut away 1 weird runs which are further than 10 sigma away from mean
    KeepIndex = (abs(KrPos-mean(KrPos))./KrPosErr)<10;
    Time = Time(KeepIndex);
    KrPos = KrPos(KeepIndex);
    KrPosErr = KrPosErr(KeepIndex);
    [par, err, chi2min, dof ] =linFit(Time,KrPos,KrPosErr);
    
    TimeMean = mean(Time);
    StartTimeStampMean = datetime('22-Oct-2019 23:17:59');
    save(savename,'Time','TimeMean','StartTimeStampMean','KrPos','KrPosErr','RunNr','par', 'err', 'chi2min', 'dof','KeepIndex');
end

qUslope = par(1);  % drift in eV/day

if strcmp(SanityPlot,'ON')
    fig1 = figure('units','normalized','pos',[0.1, 0.1,1,0.6]);
    e1 = errorbar(Time,KrPos-mean(KrPos),KrPosErr,...
        'o','Color',rgb('DodgerBlue'),'MarkerFaceColor',rgb('SkyBlue'),'LineWidth',1);
    e1.CapSize = 0;
    hold on
    pFit = plot(Time,Time.*par(1)+par(2)-mean(KrPos),'LineWidth',2,'Color',rgb('GoldenRod'));
    xlim([min(Time)-1,max(Time)+1]);
    PrettyFigureFormat('FontSize',24);
    xlabel('Time (day)')
    ylabel('Line position (eV)');
    leg = legend([e1,pFit],'Krypton (MoS)',sprintf('Fit slope %.2f \\pm %.2f meV/day , \\chi^2/dof = %.1f',par(1)*1e3,err(1)*1e3,chi2min/dof));
    leg.EdgeColor = rgb('Silver');
    leg.Location = 'southwest';
    ylim([min(KrPos-mean(KrPos)-KrPosErr)-0.02, max(KrPos-mean(KrPos)+KrPosErr)+0.02]);
end

end
