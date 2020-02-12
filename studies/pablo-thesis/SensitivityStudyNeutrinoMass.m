clear;

addpath(genpath('../../../Samak2.0'));


% Screening corrections does not exist in new TBD (commented at the end, taken from TritiumSpectrumeV)

% Object that will be used as Monte Carlo Data
A_MCD = ref_sensitivity();
A_MCD.Q_i = 18575;

% Object that will be used as Model for the Fit
Afit = ref_sensitivity();
Afit.ComputeTBDDS(); Afit.ComputeTBDIS();

% Number of Monte Carlo Simulations
Ntrials = [10000];

fitter = 'minuit';
chi2name = 'chi2Stat';


pulls = [inf,inf(1,A_MCD.nPixels*2+3)];

for trials = 1:length(Ntrials)
    
N = Ntrials(trials);
    
ResultsTable = zeros(N,10);

parfor mc = 1:N
    Data = GetData(A_MCD);
    
    F = FITC('SO',Afit,'DATA',Data,'fitter',fitter,...
        'chi2name',chi2name,...
        'pulls',pulls,...
        'fixPar','5 6',...
        'exclDataStart',1);
    % mnu, e0, bck, norm
    mnu = F.RESULTS{1}(1)+Afit.mnuSq_i;
    e0 = F.RESULTS{1}(2)+Afit.Q_i;
    bck = F.RESULTS{1}(3)+Afit.BKG_RateAllFPDSec;
    norm = F.RESULTS{1}(4)+1;
    ResultsTable(mc,:) = [mnu,e0,bck,norm,F.RESULTS{2}(1:4),F.RESULTS{3},F.RESULTS{5}];
    disp(mc)
end

save(['MC_DR',num2str(N),'.mat'],'ResultsTable')
end

function Data = GetData(A_MCD)

% Calculate Diff and Int Spectrum and add fluctuations
A_MCD.ComputeTBDDS(); A_MCD.ComputeTBDIS();

% Add the corrections
A_MCD.AddStatFluctTBDIS();

Data = [A_MCD.qU, A_MCD.TBDIS, A_MCD.TBDISE];

end
