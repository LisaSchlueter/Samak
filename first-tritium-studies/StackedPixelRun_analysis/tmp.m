addpath(genpath('../../../Samak2.0'));

% Data / Sim
mode = 'Data';
ELossFlag = 'CW_GLT'; % 'Aseev','Abdurashitov','CW_GLT'

% If simulation
StatFluct  = 'ON';
SysFluct   = 'OFF';
nSamples   = 1;

% Run List
%RunList = 'StackCD100_3hours';
RunList = 'StackCD100all';

% Range
ranges = [7]; % 9 = -200 eV; 7 = -400 eV; 1 = -1600 eV
Nranges = length(ranges);

% Choose fixed parameters (optional)
fixPar   = '1 5 6';

% Choose fitter
fitter = 'minuit';

% Choose CM
chi2name = 'chi2Stat';
%chi2name = 'chi2CMShape';


for range = 1:Nranges
    % Choose data range to analyze
    dataStart = ranges(range);
   
    switch mode
        case 'Data'
            
            MRA = MultiRunAnalysis('AnaFlag','StackPixel',...
                'chi2',chi2name,'RunList',RunList,...
                'fitter',fitter,'exclDataStart',dataStart,...
                'fixPar',fixPar,...
                'DataEffcorr','RunSummary','pullFlag',3,...
                'ELossFlag',ELossFlag);
            
            if ~strcmp(chi2name,'chi2Stat')
                MRA.ComputeCM('DataDriven','OFF','WGTS_TASR_RelErr',0.005);
            end
            
        case 'Sim'
            
            MRA = MultiRunAnalysis('AnaFlag','StackPixel',...
                'chi2',chi2name,'RunList',RunList,...
                'fitter',fitter,'exclDataStart',dataStart,...
                'fixPar',fixPar,...
                'DataEffcorr','OFF','pullFlag',3,...
                'ELossFlag',ELossFlag);
            
            if ~strcmp(chi2name,'chi2Stat')
                MRA.ComputeCM('DataDriven','OFF','WGTS_TASR_RelErr',0.005);
            end
            
            %Compute MC Data set
            switch StatFluct
                case 'ON'
                    if strcmp(SysFluct,'OFF')
                        TBDIS_Sim = mvnrnd(MRA.ModelObj.TBDIS,MRA.ModelObj.TBDIS',nSamples)'; % nSamples simulated integral spectra
                    elseif strcmp(SysFluct,'ON')
                        TBDIS_Sim = mvnrnd(MRA.ModelObj.TBDIS,MRA.FitCM+diag(MRA.ModelObj.TBDIS),nSamples)';
                    end
                case 'OFF'
                    if strcmp(SysFluct,'OFF')
                        TBDIS_Sim = repmat(MRA.ModelObj.TBDIS,1,nSamples);
                    elseif  strcmp(SysFluct,'ON')
                        TBDIS_Sim = mvnrnd(MRA.ModelObj.TBDIS,MRA.FitCM,nSamples)';
                    end
            end
            MRA.RunData.TBDIS   = TBDIS_Sim;
            MRA.RunData.TBDISE  = sqrt(TBDIS_Sim);
    end
    
    MRA.Fit();MRA.PlotFit;
    %MRA.PlotStackingModel;
end
