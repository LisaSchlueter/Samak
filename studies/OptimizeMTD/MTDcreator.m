function [qU, qUfrac, TD] = MTDcreator(varargin)

addpath(genpath('../../../Samak2.0'));

% --------------------------- PARSER START ------------------------------%
p = inputParser;
p.addParameter('E0',18575);
p.addParameter('mnuSq_i',(.5)^2,@(x)isfloat(x));
p.addParameter('Range',60);
p.addParameter('BqU',[5 10 20]); 
p.addParameter('MACE_Ba_T',6e-4);     
p.addParameter('BsBmRed',0.7); 
p.addParameter('RunTime',3*365*86400); 
p.addParameter('NuMassFactor',0.6); 
p.addParameter('BkgFraction',0.15); 
p.addParameter('Save','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('MTD_Plot','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('NuMassSignal_Plot','OFF',@(x)ismember(x,{'ON','OFF'}));
p.addParameter('AnchorBkg6G',''); % if empty, Anchor based und ROIlow is used 335=FT, [26keV, 32 keV]ROI, 238 Nikolaus PhD (70%bfields)
p.addParameter('FPD_ROIlow',14,@(x)isfloat(x));
p.addParameter('BKG_RateSec','',@(x)isfloat(x) || isempty(x)); %if empty, filled according to Ba

p.parse(varargin{:});
E0               = p.Results.E0;        % Endpoint, eV
mnuSq_i          = p.Results.mnuSq_i;   % nu mass squared, eV^2
Range            = p.Results.Range;     % depth, eV
BqU              = p.Results.BqU;       % background
MACE_Ba_T        = p.Results.MACE_Ba_T; % Ba in T
BsBmRed          = p.Results.BsBmRed;   % WGTS/Max B field reduction 
RunTime          = p.Results.RunTime;   % Run time, seconds
NuMassFactor     = p.Results.NuMassFactor;    % Run time, seconds
BkgFraction      = p.Results.BkgFraction;     % Run time, seconds
Save             = p.Results.Save;            % 
MTD_Plot         = p.Results.MTD_Plot;         % 
NuMassSignal_Plot= p.Results.NuMassSignal_Plot;% 
AnchorBkg6G      = p.Results.AnchorBkg6G;
FPD_ROIlow       = p.Results.FPD_ROIlow;
BKG_RateSec      = p.Results.BKG_RateSec;
% --------------------------- PARSER ENDS  ------------------------------%

% ---------------------- Simulation Model  ------------------------------%
WGTS_B_T    = BsBmRed*3.6;
MACE_Bmax_T = (WGTS_B_T/3.6)*6;

if isempty(BKG_RateSec)
BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,'FPD_ROIlow',FPD_ROIlow,...
    'AnchorBkg6G',AnchorBkg6G);
end
%BKG_RateSec = GetBackground('MACE_Ba_T',MACE_Ba_T,'WGTS_B_T',WGTS_B_T,'Anchor6G',Anchor6G,'Anchor6GValue',Anchor6GValue);
%BKG_RateSec = 10e-3;
% -----------------------------------------------------------------------%


% -------------------- MTD Step 1: Build Flat MTD -----------------------%
[SqU , SqUfrac ]= FlatMTD_S_creator('Range',Range,'Save','OFF');
qU = SqU; qUfrac=SqUfrac ; 
TD = 'MTDcreatorTMP'; save(TD,'TD','qU','qUfrac');
if ismac
    !mv MTDcreatorTMP.mat ../../simulation/katrinsetup/TD_DataBank/
else
    system('mv MTDcreatorTMP.mat ../../simulation/katrinsetup/TD_DataBank/');
end
comopt = {'BKG_RateAllFPDSec',BKG_RateSec,'MACE_Ba_T',MACE_Ba_T,...
    'WGTS_B_T',WGTS_B_T,'MACE_Bmax_T',MACE_Bmax_T,...
    'TimeSec',RunTime};

% -------------------- MTD Step 2: Build IsoStat MTD --------------------%
S = ref_TBD_NominalKATRIN(comopt{:},'TD',TD,'mnuSq_i',0);
S.ComputeTBDDS; S.ComputeTBDIS; SqUfrac = qUfrac ./ (S.TBDIS); SqUfrac = SqUfrac./sum(SqUfrac);
SqUfrac  = (S.qU<(E0-30)).*(SqUfrac .* 5) + (S.qU>=(E0-30)).*(SqUfrac .* 1);
SqUfrac  = SqUfrac./sum(SqUfrac); 
S.qUfrac  = SqUfrac; 
S.ComputeTBDDS; S.ComputeTBDIS;

% -------------------- MTD Step 3: Build IsoStat+ MTD -------------------%
S1 = ref_TBD_NominalKATRIN(comopt{:},'TD',TD,'mnuSq_i',mnuSq_i);
S1.qUfrac  = SqUfrac; S1.ComputeTBDDS; S1.ComputeTBDIS;
% signal  = (S1.TBDIS-S.TBDIS).^2./S1.TBDIS.*(S.qU>=(E0-30) & S.qU<=(E0-1));
s=S.TBDIS-(S.qUfrac.*S.TimeSec.*S.BKG_RateSec);
b=S.qUfrac.*S.TimeSec.*S.BKG_RateSec;
Esb2=(s.^2+b.^2);
signal  = (S1.TBDIS-S.TBDIS).^2./Esb2.*(S.qU>=(E0-30) & S.qU<=(E0-4));
%signal  = (S.TBDIS-S1.TBDIS)./sqrt(S.TBDIS).*(S.qU>=(E0-30));
SqUfracNuMass = signal./sum(signal);
SqUfrac = (1-NuMassFactor) .* SqUfrac + (NuMassFactor) .* SqUfracNuMass;
SqUfrac=SqUfrac./sum(SqUfrac); 
% -----------------------------------------------------------------------%
% disp(sum(SqUfrac))
SqUfrac=(SqUfrac<=1e-3).*1e-3 + SqUfrac.*(SqUfrac>1e-3);
SqUfrac=SqUfrac./sum(SqUfrac); 
% disp(sum(SqUfrac));
% -----------------------------------------------------------------------%
S.qUfrac= SqUfrac; S1.qUfrac = SqUfrac;  
% -----------------------------------------------------------------------%

% ------------------ MTD Step 4: Flat Background  -----------------------%
[BqU BqUfrac ]= FlatMTD_B_creator('BqU',BqU,'Save','OFF');
qU     = [SqU ; BqU'];
qUfrac = 1 ./ numel(qU);
TD = 'MTDcreatorTMP';
save(TD,'TD','qU','qUfrac');
if ismac
!mv MTDcreatorTMP.mat ../../simulation/katrinsetup/TD_DataBank/
else
    system('mv MTDcreatorTMP.mat ../../simulation/katrinsetup/TD_DataBank/');
end
% -----------------------------------------------------------------------%

% --------- MTD Step 5: IsoStat + NuMass Background  --------------------%

% Renormalize IsoStat + NuMass MTD
SqUfrac = (1-BkgFraction).*SqUfrac;
% Renormalize Background MTD
BqUfrac = BkgFraction.*BqUfrac;
% Add IsoStat + NuMass MTD + Background MTD
qU = [SqU;BqU'];
qUfrac = [SqUfrac ; BqUfrac];
TD = 'MTDcreator';
% Save
save([TD '.mat'],'TD','qU','qUfrac');
save([TD '.mat'],'TD','qU','qUfrac');
% Save .txt File
qUsamak=qU';
qUfracsamak=qUfrac';
save([TD '.txt'],'qUsamak','qUfracsamak','-ascii')

% Save in DataBank
filename = sprintf('MTDcreator_E0%.1f_%.0feV_B%0.f_Ba%.1f_RedF%0.1f_NuMF%0.2f_BkgF%0.2f_B%.0f.mat',...
            E0,Range,sum((BqU-E0)),MACE_Ba_T/1e-4,BsBmRed,...
            NuMassFactor,BkgFraction,round(BKG_RateSec*1000));
system(['mv MTDcreator.mat ../../simulation/katrinsetup/TD_DataBank/' filename]);

% Plot
switch MTD_Plot
    case 'ON'
        PlotMTD('TDprefix','TDCreator','RunTimeSec',RunTime)
end

% Plot Spectrum / Signal
switch NuMassSignal_Plot
    case 'ON'
        figure(2)
        subplot(2,1,1)
        S  = ref_TBD_NominalKATRIN(comopt{:},'TD',TD,'mnuSq_i',0);
        S.ComputeTBDDS; S.ComputeTBDIS;
        S1 = ref_TBD_NominalKATRIN(comopt{:},'TD',TD,'mnuSq_i',mnuSq_i);
        S1.ComputeTBDDS; S1.ComputeTBDIS;
        errorbar(S.qU-S.Q,(S1.TBDIS./S.TBDIS),S1.TBDISE./S.TBDIS,'LineWidth',1,'Color','Black','LineStyle','-');
        ylabel('Ratio','FontSize',14);
        xstr = sprintf('E-%.1f (eV)',E0);
        xlabel(xstr,'FontSize',14);
        PrettyFigureFormat
        subplot(2,1,2)
        plot(S.qU-S.Q,S1.TBDISE./S.TBDIS,'LineWidth',1,'Color','Black','LineStyle','-');
        xlabel(xstr,'FontSize',14);
        ylabel('Relative Error','FontSize',14);
        PrettyFigureFormat
end

% Remove Temporary Files
if ismac
    !rm ../../simulation/katrinsetup/TD_DataBank/MTDcreatorTMP.mat
else
    system('rm ../../simulation/katrinsetup/TD_DataBank/MTDcreatorTMP.mat')
end
end

