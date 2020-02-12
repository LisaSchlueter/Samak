function [PnoEE, fit_p_data_ee, fit_p_data_noee]=RunWiseFSDeeonsetRejection(varargin)
%  Compute the probability to reject the no electronics excited states
%  hypothesis based on a frequentist set of 'nfits' simulations
%  n Fits of First Tritium stacked-pixel single run for
%  - E0 , m , B , N
%  - Step 1) Pgs is fixed to its expected value
%  - Step 2) Pgs is fixed to 1 (noEE)
%
% Input:
%  - run             run number (first tritium)
%  - nfit            number of fits
%  - runSim          ON=perform full simulation /OFF=read old file
%  - exclDataStart:  7=400eV, 9=200eV, ... 
% Output:
% - PnoEE            probability for no EE
% - fit_p_data_ee    fit parameters - Pg fixed to nominal
% - fit_p_data_noee  fit parameters - Pg fixed to 1
% 
% T. Lasserre, Pablo Morales, Lisa Schlueter
% Last Update: June 24 2018
%

% Parser
p = inputParser;
p.addParameter('run',40680);
p.addParameter('nfit',50,@(x)isfloat(x));
p.addParameter('exclDataStart',7,@(x)isfloat(x));
p.addParameter('runSim','ON',@(x)ismember(x,{'ON','OFF'}));

p.parse(varargin{:});

run            = p.Results.run;
nfit           = p.Results.nfit;
exclDataStart  = p.Results.exclDataStart;
runSim         = p.Results.runSim;

tmpstr  = sprintf('fsdonset-run%gexcl%g.mat',run,exclDataStart);
switch runSim
    case 'OFF'
%    tmpstr  = sprintf('fsdonset-run40683excl%g.mat',exclDataStart);
    tmpstr  = sprintf('fsdonset-run40683excl7.mat');
    tmpfile = load(tmpstr);
    fit_p_sim_ee   = tmpfile.fit_p_sim_ee;
    fit_p_sim_noee = tmpfile.fit_p_sim_noee;
    case 'ON'
        % Simulation - FSD - with EE
        fit_p_sim_ee    =  FT_FitFSD('run',run,'FSD','eeON','Mode','Sim','StatFluc','ON','nfit',nfit,'displayFit','OFF','exclDataStart',exclDataStart);
        % Simulation FSD   - No EE
        fit_p_sim_noee  =  FT_FitFSD('run',run,'FSD','eeOFF','Mode','Sim','StatFluc','ON','nfit',nfit,'displayFit','OFF','exclDataStart',exclDataStart);
        % Save 
        save(tmpstr,'fit_p_sim_ee','fit_p_sim_noee');
end

%% Data -  - FSD - with EE
fit_p_data_ee   =  FT_FitFSD('run',run,'FSD','eeON','Mode','Data','nfit',1,'exclDataStart',exclDataStart);
%% Data -  - FSD - No EE
fit_p_data_noee =  FT_FitFSD('run',run,'FSD','eeOFF','Mode','Data','nfit',1,'exclDataStart',exclDataStart);

%% No-EE Hypothesis Rejection
PnoEE=100-sum(fit_p_sim_ee(:,13)>fit_p_data_noee(13))/numel(fit_p_sim_ee(:,13))*100;

%% Plot Simulation - FSD - with EE
fig1 = figure('Name','VFT Fits','NumberTitle','off','rend','painters','pos',[10 10 800 600]);
hee=histogram(fit_p_sim_ee(:,13),'Facecolor',rgb('CadetBlue'),'FaceAlpha',.5,'EdgeColor','none');
xlabel('\chi^2');
ylabel('simulations / data');
xlim([0 120]);
hold on
hnoee=histogram(fit_p_sim_noee(:,13),'Facecolor',rgb('Khaki'),'FaceAlpha',.5,'EdgeColor','none');
lnoee=line([fit_p_data_noee(13) fit_p_data_noee(13)],[0 max(hnoee.Values)],'LineWidth',2,'Color',rgb('Khaki'),'LineStyle','--');
lee=line([fit_p_data_ee(13) fit_p_data_ee(13)],[0 max(hee.Values)],'LineWidth',2,'Color',rgb('CadetBlue'),'LineStyle','--');
hold off
titlestr=sprintf('run %g Simulation / Data - Electronic Excitations (EE) Onset',run);
title(titlestr)
box off
lgdstr1=sprintf('%g simulations: fit with EE',nfit);
lgdstr2=sprintf('%g simulations: fit without EE',nfit);
lgdstr4=sprintf('%g simulations: fit without EE \n Hypothesis rejected at %g %% CL',nfit,PnoEE);
lgd=legend([hee,hnoee,lee,lnoee],...
    {lgdstr1,lgdstr2,'Data: fit with EE',lgdstr4});
lgd.FontSize = 12;
legend boxoff
PrettyFigureFormat
strfile=sprintf('figures/fsdontset_run%g.pdf',run);
publish_figurePDF(fig1,strfile);

