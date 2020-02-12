%
% Samak KATRIN Simulation and Analysis
% 
% Fit Run List
%  - Stacked-Pixels
%  - Stacked-Runs
%
% Input:
% - RunList:      0 for all, or a vector like [run1 run2 ...]
% - chi2:         chi2CM, chi2Stat,...
% - Mode          Sim, Data
% Output:
% - Fit and Plot
% 
% T. Lasserre, Pablo Morales, Lisa Schlueter
% Last Update: June 24 2018
%

close all

% Option: Run List
RunList = 0; 

% Option: chisquare
%chi2     = 'chi2Stat';
chi2    = 'chi2CMShape';
%chi2    = 'chi2CM';

% Option: data/sim
Mode     = 'Data';
%Mode     = 'Sim';

% Fix Run List
if RunList == 0
    [ n, datafiles , RunList ] = GetRunList( '../../tritium-data/mat/','*ex2.mat',1,'string');
%       RunList1=RunList(RunList>40538 & RunList<40666);
%       RunList2=RunList(RunList>40666 & RunList<40694);
%       RunList= [ RunList1 RunList2]; 
RunList=RunList(RunList>40666 & RunList<40685);
RunList=RunList(RunList>40666 & RunList<40694);
RunList = 'StackCD100';
end

% MultiRunAnalysis
FT = MultiRunAnalysis('RunList',RunList,'chi2',chi2,'exclDataStart',4,...
    'fixPar','1 5 6','ringCutFlag','ex2','AnaFlag','StackPixel','fitter','matlab',...
    'DataEffCorr','OFF');

% Option for simulation
switch Mode
    case 'Sim'
        FT.RunData.TBDIS = zeros(FT.ModelObj.nqU,1);
        FT.LoadSingleRunObj();
        for i=1:1:numel(RunList)
            %FT.SingleRunObj{i}.AddStatFluctTBDIS;
            FT.RunData.TBDIS  = FT.RunData.TBDIS + FT.SingleRunObj{i}.TBDIS;
        end
        FT.RunData.qU     = FT.ModelObj.qU;
end

% chisquare option - Compute Covariance Matrix
if ~strcmp(chi2,'chi2Stat')
    FT.ComputeCM('DataDriven','OFF','WGTS_CD_MolPerCm2_RelErr',0.05);
    %FT.ComputeCM();
    FT.ComputeCM;
end

% Fit / Plot
FT.Fit('CATS','OFF');
FT.PlotFit('saveplot','ON');
FT.PlotDataModel


