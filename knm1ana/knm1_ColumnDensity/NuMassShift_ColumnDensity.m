% Study of the impact of a column density bias to the neutrino mass 
% squared central value
chi2 = 'chi2CMShape';
Mode = 'Twin';
RhoDbias = 0.01;

switch Mode
    case 'Real'
        A = MultiRunAnalysis('RunList','KNM1','chi2','chi2Stat','DataType','Real','minuitOpt','min  ; migrad ');
        A.fixPar = '1 5 6 7 8 9 10 11';
    case 'Twin'
        A = MultiRunAnalysis('RunList','KNM1','chi2','chi2Stat','DataType','Twin','minuitOpt','min  ; migrad ');
        A.fixPar = '5 6 7 8 9 10 11';
    case 'Asimov'
        A = MultiRunAnalysis('RunList','KNM1','chi2','chi2Stat','minuitOpt','min  ; migrad ; imp;  minos;');
        A.fixPar = '5 6 7 8 9 10 11';
        % Get Background and Normalization from Data
        A.InitModelObj_Norm_BKG;
        Bkg = A.ModelObj.BKG_RateSec;
        
        % Calculate Asimov Data
        A.ModelObj.ComputeTBDDS;
        A.ModelObj.BKG_RateSec = Bkg;
        A.ModelObj.ComputeTBDIS;
        
        % Fill Data with Asimov
        A.RunData.TBDIS = A.ModelObj.TBDIS;
        A.RunData.TBDISE = sqrt(A.ModelObj.TBDIS);
end

CD_i = A.ModelObj.WGTS_CD_MolPerCm2;
%% Fit

A.exclDataStart = 2;
E0 = zeros(3,1);
N = zeros(3,1);

A.SimulateStackRuns('WGTS_CD_MolPerCm2',CD_i*(1-RhoDbias));
A.Fit
mNuSq(1) = A.FitResult.par(1);
E0(1) = A.FitResult.par(2);
N(1) = A.FitResult.par(4);

A.SimulateStackRuns('WGTS_CD_MolPerCm2',CD_i);%,'BKG_RateAllFPDSec',BKG);
A.Fit; A.PlotFit;
mNuSq(2) = A.FitResult.par(1);
E0(2) = A.FitResult.par(2);
N(2) = A.FitResult.par(4);

A.SimulateStackRuns('WGTS_CD_MolPerCm2',CD_i*(1+RhoDbias));%,'BKG_RateAllFPDSec',BKG);
A.Fit; A.PlotFit;
mNuSq(3) = A.FitResult.par(1);
E0(3) = A.FitResult.par(2);
N(3) = A.FitResult.par(4);

% fit with different column densities
DataType = 'Twin';
RunList = 'KNM1';
chi2name = 'chi2Stat';

exclDataStart = 14;
CDRel = [0.98,0.99,1,1.01,1.02]; %relative column densities
nCD   = numel(CDRel);
savedir = [getenv('SamakPath'),'knm1ana/knm1_ColumnDensity/results'];
savename = [savedir,sprintf('NuMassShift_ColumnDensity_FitResults_%s_%s_excl%.0f_%.0f-%.0f_%.0fnCD_%s.mat',...
    DataType,RunList,exclDataStart,min(CDRel)*100,max(CDRel)*100,nCD,chi2name)];

if ~exist(savedir,'dir')
    system(['mkdir -p ',savedir]);
end

if exist(savename,'file')
    load(savename);
else
    A = MultiRunAnalysis('RunList',RunList,'chi2','chi2Stat','DataType',DataType);
    A.chi2 = chi2name;
    A.fixPar = '5 6 7 8 9 10 11';
    A.exclDataStart = exclDataStart;
    CD_i = A.ModelObj.WGTS_CD_MolPerCm2;
    
    %% Fit
    mNuSq = zeros(nCD,1);
    E0    = zeros(nCD,1);
    N     = zeros(nCD,1);
    
    for i=1:nCD
        A.SimulateStackRuns('WGTS_CD_MolPerCm2',CD_i*CDRel(i));
        A.Fit
        mNuSq(i) = A.FitResult.par(1);
        E0(i) = A.FitResult.par(2);
        N(i) = A.FitResult.par(4);
    end
end
fprintf('---- nu-mass shifts --- \n');
fprintf(' %.2g \n',mNuSq);

%%
 t = PrintTable('Neutrino Mass Shift Induced by Column Density Bias');
            t.addRow('Column Density','Neutrino Mass Squared Shift (eV$^2$)','Endpoint Shift (eV)');
            t.addRow(sprintf('Nominal - %.2f percent',RhoDbias*100),round(mNuSq(1),2),round(E0(1),3));
            t.addRow('Nominal',round(mNuSq(2),2),round(E0(2),3));
            t.addRow(sprintf('Nominal + %.2f percent',RhoDbias*100),round(mNuSq(3),2),round(E0(3),3));           
            t.addRow('','','');
            t.display;
            %t.HasHeader = true;
            t.Format = 'tex';
            if strcmp(chi2,'chi2Stat')
            t.Caption = sprintf('Neutrino Mass Shift Induced by Column Density Bias of %.1f percent - %.0f runs (stat)',RhoDbias*100,numel(A.RunList)');
            else
            t.Caption = sprintf('Neutrino Mass Shift Induced by Column Density Bias of %.1f percent - %.0f runs (stat+sys)',RhoDbias*100,numel(A.RunList)');
            end                
            t.print;