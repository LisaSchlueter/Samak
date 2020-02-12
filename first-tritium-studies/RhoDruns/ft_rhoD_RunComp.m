%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intro info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Runs with 48.08% CD = [40763,40764,40765,40766]
% Runs with 24.41% CD = [40794:40805]
% Runs with 71.94% CD = [40926:40935]
% Runs with 100% CD   = [40667,40668,40669,40670,40671,40672,40673,
%                        40674,40675,40676,40677,40678,40679,40680,
%                        40681,40682,40683,40684,40685,40686,40687,...
%                        40688,40689,40690,40691,40692,40693];

clear;

ComparisonPlots = 'OFF';

% Arrange runs in cells
RCD{4} = [40668]; % 100
RCD{2} = [40763,40764,40765,40766]; % 48
RCD{1} = [40794:40804]; %24
RCD{3} = [40926:40935]; %71

% Allocate memory
MRA = cell(4,1);
Stacked = cell(4,1);
CD = zeros(4,1);
CDall = zeros(15,4);

% Check the CD of each run individually, to confirm that the groups of runs
% were chosen correctly
for rr = 1:4
    Runs = RCD{rr};
    for r = 1:length(Runs)
        R = load([num2str(Runs(r)),'ex2.mat']);
        CDall(r,rr) = R.WGTS_CD_MolPerCm2;
    end
end

% Loop to extract data
TTc = zeros(4,1);
for r = 1:4
    MRA{r} = MultiRunAnalysis('RunList',RCD{r},'chi2','chi2Stat','exclDataStart',9,...
    'fixPar','1 5 6','ringCutFlag','ex2','AnaFlag','StackPixel');
    Stacked{r} = MRA{r}.StackedRuns;
    CD(r) = MRA{r}.RunData.WGTS_CD_MolPerCm2;
    TimeSecFrac(:,r) = MRA{r}.RunData.TimeSec.*MRA{r}.RunData.qUfrac;
    TBDIS(:,r) = MRA{r}.RunData.TBDIS;
    qU(:,r) = MRA{r}.RunData.qU;
    DTc(r) = MRA{r}.RunData.WGTS_MolFrac_DT;
    HTc(r) = MRA{r}.RunData.WGTS_MolFrac_HT;
    %TTc(r) = MRA{r}.RunData.WGTS_MolFrac_TT;
    Tp(r) = 0.5*DTc(r) + 0.5*HTc(r) + TTc(r);
    
end

if strcmp(ComparisonPlots,'ON')

% Normalize to one of the runs (in this case to run 40668, with 100% CD)
CDRatio = CD/CD(4);
TimeRatio = TimeSecFrac./TimeSecFrac(:,4);
TritiumPurityRatio = Tp./Tp(4);

%% Plot taking into account only the normalization with CD and time
figure(1)
for r = 1:4
    hold on
    %if r == 1; TBDIS(:,r) = 2*TBDIS(:,r); end
    TBDISreduced(:,r) = (TBDIS(:,r)./TimeRatio(:,r))./CDRatio(r)./1;
    plotCD = plot(qU(:,r)-18575,TBDISreduced(:,r));
    plotCD.LineWidth = 2;
    set(gca, 'YScale', 'log')
    legcell{r} = ['\rho d =  ',num2str(round(CD(r),3,'significant')),'cm^{-2} (',...
        num2str(round(CDRatio(r)*1e2,3,'significant')),'%)']; 
    hold off
end

title('Integral Spectra normalized to Run 40668 (Time and \rho d)');
ylabel('Normalized Counts'); 
xlabel('E-E_{0 (18575)} (eV)')
axis([-1600 50 0 1e8 ]);
PrettyFigureFormat;
legend(legcell,'Location','Best')
export_fig(['plots/','Comp1.pdf'],'-pdf');
close all

% Plot taking into account only the normalization with CD, time, and tritium purity
figure(2)
for r = 1:4
    hold on
    TBDISreduced(:,r) = (TBDIS(:,r)./TimeRatio(:,r))./CDRatio(r)./TritiumPurityRatio(r);
    plotCD = plot(qU(:,r)-18575,TBDISreduced(:,r));
    plotCD.LineWidth = 2;
    set(gca, 'YScale', 'log')
    legcell{r} = ['\rho d =  ',num2str(round(CD(r),3,'significant')),'cm^{-2} (',...
        num2str(round(CDRatio(r)*1e2,3,'significant')),'%)']; 
    hold off
end

title('Integral Spectra normalized to Run 40668 (Time, \rho d and Tritium purity)');
ylabel('Normalized Counts');
xlabel('E-E_{0 (18575)} (eV)')
axis([-1600 50 0 1e8]);
PrettyFigureFormat;
legend(legcell,'Location','Best');
export_fig(['plots/','Comp2.pdf'],'-pdf');


%% RATIO - Plot taking into account only the normalization with CD, time, and tritium purity
figure(3)
for r = 1:4
    hold on
    TBDISreduced(:,r) = (TBDIS(:,r)./TimeRatio(:,r))./CDRatio(r)./TritiumPurityRatio(r);
    plotCD = plot(qU(:,r)-18575,TBDISreduced(:,r)./TBDISreduced(:,4));
    plotCD.LineWidth = 2;
    set(gca, 'YScale', 'lin')
    legcell{r} = ['\rho d =  ',num2str(round(CD(r),3,'significant')),'cm^{-2} (',...
        num2str(round(CDRatio(r)*1e2,3,'significant')),'%)']; 
    hold off
end

title('Integral Spectra normalized to Run 40668 (Time, \rho d and Tritium purity)');
ylabel('Ratio');
xlabel('E-E_{0 (18575)} (eV)')
xlim([-1600 50]);
PrettyFigureFormat;
legend(legcell,'Location','Best');
export_fig(['plots/','Comp3.pdf'],'-pdf');



if ismac; unix('touch plots/SummaryRhoDHandle.pdf'); end
 
FileNames = {['plots/','Comp1.pdf'],['plots/','Comp2.pdf'],['plots/','Comp3.pdf']};

append_pdfs(['plots/SummaryRhoDHandle.pdf'],FileNames{:});
end




