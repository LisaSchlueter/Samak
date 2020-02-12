function [] = comparespectra(varargin)
%%
%
% P. I. Morales Guzmán, December 2017
%

%% Constants
c   = 299792458;    % speed of light
LW  = 1;            % linewidth for plots
%% Parser
p = inputParser;

% File names from the curves to be compared (Data from SAMAK can be a
% 2 or 3-column (qU rate uncertainty, for the integral spectrum) matrix or a TBD object)
p.addParameter('SAMAKdata','',@(x) isa(x,'TBD') || isa(x,'double'));
p.addParameter('SSCdata','',@(x) isa(x,'double'));

% Type of curve to compare (phase space, differential spectrum, integral
% spectrum, Fermi function, response function)
p.addParameter('compType','DIFF',@(x)ismember(x,{'PS','DIFF','INT','Fermi','RF'}));
p.addParameter('PhaseSpaceType','direct',@(x)ismember(x,{'direct','elementbyelement'}));

% Normalization flag (apply normalization to spectra)
p.addParameter('normFlag','ON',@(x)ismember(x,{'ON','OFF'}));

% Plot flag and text options
p.addParameter('plot_title','',@(x) isa(x,'char'));
p.addParameter('plot_save','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('plot_file','.png',@(x)ismember(x,{'.png','.pdf','.eps','.fig'}));
p.addParameter('spec1','SAMAK',@(x) isa(x,'char'));
p.addParameter('spec2','SSC',@(x) isa(x,'char'));
p.addParameter('showinfo','ON',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('errorbarFlag','ON',@(x)ismember(x,{'ON','OFF'}));



% Parser
p.parse(varargin{:});
SAMAKdata       = p.Results.SAMAKdata;
SSCdata         = p.Results.SSCdata;
compType        = p.Results.compType;
PhaseSpaceType  = p.Results.PhaseSpaceType;
normFlag        = p.Results.normFlag;
plot_title      = p.Results.plot_title;
plot_save       = p.Results.plot_save;
plot_file       = p.Results.plot_file;
spec1           = p.Results.spec1;
spec2           = p.Results.spec2;
showinfo        = p.Results.showinfo;
errorbarFlag    = p.Results.errorbarFlag;



%% Import data

% If the comparison is for the phase space and element by element, the
% spectra to compare are given later, not directly from a file
if ~(strcmp(compType,'PS') && strcmp(PhaseSpaceType,'elementbyelement'))
    TBD_SSC = SSCdata;
end

if isa(SAMAKdata,'TBD')
    switch compType
        case 'PS'
            TBD_SAMAK(:,1) = SAMAKdata.Te;
            switch PhaseSpaceType
                case 'direct'
                    TBD_SAMAK(:,2) = SAMAKdata.PhaseSpace;
                case 'elementbyelement'
                    TBD_SAMAK(:,2) = readPhaseSpaceElementByElement(SAMAKdata);
            end
        case 'DIFF'
            TBD_SAMAK(:,2) = SAMAKdata.TBDDS;
            TBD_SAMAK(:,1) = SAMAKdata.Te;
        case 'INT'
            TBD_SAMAK(:,2) = SAMAKdata.TBDIS;
            TBD_SAMAK(:,1) = SAMAKdata.qU;
            TBD_SAMAK(:,3) = SAMAKdata.TBDISE;
        case 'Fermi'
            TBD_SAMAK(:,2) = SAMAKdata.ComputeFermiCorr();
            TBD_SAMAK(:,1) = SAMAKdata.Te;
        case 'RF'
            TBD_SAMAK(:,2) = SAMAKdata.KTF(SAMAKdata.Te(min(TBD_SSC(:,1),max(TBD_SSC(:,1)))),...
                18575,SAMAKdata.MACE_R_eV);
            TBD_SAMAK(:,1) = SAMAKdata.Te;
    end
else
    TBD_SAMAK = SAMAKdata;
end

% For the phase space from SSC, just the values of the function are given
% and the range sould be copied from SAMAK.
if strcmp(compType,'PS') && strcmp(PhaseSpaceType,'elementbyelement')
    TBD_SSC = readPhaseSpaceElementByElement(1);
    TBD_SSC(:,1) = TBD_SAMAK(:,1);
end

%% Calculate normalization factor
    NormFactSAMAK = 1;
    NormFactSSC = 1;
if strcmp(normFlag,'ON') && ~(ismember(compType,{'Fermi','RF'}))
    BckgSAMAK   = 0;
    BckgSSC     = 0;
    switch compType
        case 'PS'
            if TBD_SAMAK(1,2) > 1e6*TBD_SSC(1,2)
                NormFactSAMAK = 1/c;
            else
                NormFactSAMAK = 1;
            end
            NormFactSSC   = 1;
        case 'DIFF'
            NormFactSAMAK = 1/c;
            NormFactSSC   = 3.276645539579804e9/0.004557704049428*...
                (TBD_SAMAK(1,2)/TBD_SSC(1,2)); %A.NormFactorTBDDS
        case 'INT'
            NormFactSAMAK = 1;
            NormFactSSC   = TBD_SAMAK(1,2)/TBD_SSC(1,2);
            BckgSAMAK     = TBD_SAMAK(end,2);
            BckgSSC       = TBD_SSC(end,2);
    end
    % Normalization should only be applied to the spectrum without
    % background, since normally a flat background is used, and it is the
    % same regardless of differences on the spectrum.
    TBD_SSC(:,2)     = (TBD_SSC(:,2) - BckgSSC)*NormFactSSC + BckgSSC;
    TBD_SAMAK(:,2)   = (TBD_SAMAK(:,2) - BckgSAMAK)*NormFactSAMAK + BckgSAMAK;
end

%% Set range
% The plots won't work if the qU values are different.
if TBD_SSC(end,1) > 18575
    TBD_SSC(:,1) =  TBD_SSC(:,1) - 18575;
end

if TBD_SAMAK(end,1) > 18575
    TBD_SAMAK(:,1) =  TBD_SAMAK(:,1) - 18575;
end

%% Plot direct comparison
figure(1)
hold on
switch compType
    case {'PS','DIFF','Fermi','RF'}
        plot(TBD_SSC(:,1),TBD_SSC(:,2),'LineWidth',LW)
        plot(TBD_SAMAK(:,1),TBD_SAMAK(:,2),'LineWidth',LW)
    case 'INT'
        if strcmp(errorbarFlag,'ON')
            errorbar(TBD_SSC(:,1),TBD_SSC(:,2),TBD_SSC(:,3),'LineWidth',LW)
            errorbar(TBD_SAMAK(:,1),TBD_SAMAK(:,2),TBD_SAMAK(:,3),'LineWidth',LW)
        else
            plot(TBD_SSC(:,1),TBD_SSC(:,2),'LineWidth',LW)
            plot(TBD_SAMAK(:,1),TBD_SAMAK(:,2),'LineWidth',LW)
        end
end
hold off
PrettyFigureFormat;
legend({spec1,spec2})
title([compType,' Comparison ',plot_title])
xlabel('E - E_0 [eV]');
if ismember(compType,{'INT','DIFF'})
    ylabel('cps per eV');
elseif ismember(compType,{'RF','Fermi'})
    ylabel('probability');
else
    ylabel('a. u.')
end
annotation('textbox',[0.2 0.5 0.3 0.3],'String',...
    ['Norm. SSC = ',num2str(round(NormFactSSC,4))],'FitBoxToText','on');

if strcmp(plot_save,'ON')
    file_num = 0;
    while exist(['DC_CompFigs/BothTBDDS',compType,num2str(file_num),plot_file], 'file') == 2
        file_num = file_num + 1;
    end
    saveas(gcf,['DC_CompFigs/BothTBDDS',compType,num2str(file_num),plot_file])
end

%% Calculate differences

difference = TBD_SAMAK(:,2) - TBD_SSC(:,2);
if strcmp(compType,'INT')
    difference(abs(difference) < 5e-9) = 0;
end
average = 0.5*(TBD_SAMAK(:,2) + TBD_SSC(:,2));
relativeDifference = difference./average;
relativeDifference(isnan(relativeDifference)) = 0;

figure(2)
switch compType
    case {'PS','DIFF','Fermi','RF'}
        plot(TBD_SAMAK(:,1),relativeDifference,'LineWidth',LW);
    case 'INT'
        if strcmp(errorbarFlag,'ON')
            errorINT2 = (average).^(-4).*((TBD_SSC(:,2).*TBD_SAMAK(:,3)).^2 +...
                (TBD_SAMAK(:,2).*TBD_SSC(:,3)).^2);
            errorINT = sqrt(errorINT2);
            errorbar(TBD_SAMAK(:,1),relativeDifference,errorINT,'LineWidth',LW);
        else
            plot(TBD_SAMAK(:,1),relativeDifference,'LineWidth',LW);
        end
end
% axis([-inf 18574.7 0 0.4]);
if ismember(compType,{'INT','DIFF'})
    title(['Tritium Beta Decay ',compType,' Relative Diff. ',plot_title]);
else
    title([compType,' Relative Diff. ',plot_title]);
end

xlabel('E - E_0 [eV]');
legend([spec1,' and ',spec2]);
PrettyFigureFormat;

if strcmp(plot_save,'ON')
    file_num = 0;
    while exist(['samak/studies/datachallenge/DC_CompFigs/RelativeDifference',compType,num2str(file_num),plot_file], 'file') == 2
        file_num = file_num + 1;
    end
    saveas(gcf,['samak/studies/datachallenge/DC_CompFigs/RelativeDifference',compType,num2str(file_num),plot_file])
end

relativeDifference = relativeDifference(~isnan(relativeDifference));

if strcmp(showinfo,'ON')
    fprintf('===============================================\n');
    fprintf('Max diff. (absolute value)       = %e\n', max(abs(difference)));
    fprintf('Mean diff. (absolute value)      = %e\n', mean(abs(difference)));
    fprintf('Max rel. diff. (absolute value)  = %e\n', max(abs(relativeDifference)));
    fprintf('Mean rel. diff. (absolute value) = %e\n', mean(abs(relativeDifference)));
    fprintf('Sum of squares                   = %e\n', sum(difference.^2));
    fprintf('Normalized by hand  = %f \n', NormFactSSC);
    fprintf('===============================================\n');
end
end

function PhaseSpace = readPhaseSpaceElementByElement(varargin)

if isa(varargin{1},'TBD')
    A = varargin{1};
    PhaseSpace = [A.Te, A.pe.*(A.Te+A.me).*(A.Q-A.Te).^2.*(A.Te<18575)];
else
    EmeSSC  = importdata(['DC_Data/SSC_e_m.txt']);
    SSCp    = importdata(['DC_Data/SSCp.txt']);
    FSDSum  = importdata(['DC_Data/SSCFSDSum.txt']);
    FermiSSC= importdata(['DC_Data/SSCFermi.txt']);
    PhaseSpace(:,2) = SSCp(:,2).*EmeSSC(:,2).*FSDSum(:,2).*FermiSSC;
end

end