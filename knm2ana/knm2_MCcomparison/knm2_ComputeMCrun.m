% KNM2 Monte Carlo data challenge
% Generate 1 MC run
% Convert to HDF5-format for other fitters

TwinBias_mnuSq = 0;%0.28;
range = 90;
KTFFlag = 'WGTSMACE';%_NIS1';
fitter = 'minuit';
RecomputeFlag = 'ON';
PlotChi2Curve = 'OFF';

if ~strcmp(KTFFlag,'WGTSMACE')
    scatStr = ['_',KTFFlag];
else
    scatStr = '';
end

% label
savedir = [getenv('SamakPath'),'knm2ana/knm2_MCcomparison/results/'];
MakeDir(savedir);
savename = [savedir,sprintf('knm2_ComputeMCrunFitResult_mNuSq%.2feV2_%.0feV_%s%s.mat',TwinBias_mnuSq,range,fitter,scatStr)];

if exist(savename,'file') && strcmp(RecomputeFlag,'OFF')
    load(savename);
else
    nStack = 100;
    RunAnaArg = {'RunNr',56341,...         %
        'fixPar','mNu E0 Bkg Norm',...     % free Parameter !!
        'DataType','Twin',...              % Real, Twin or Fake
        'FSDFlag','BlindingKNM2',...       % final state distribution (theoretical calculation)
        'ELossFlag','KatrinT2',...         % energy loss function     ( different parametrizations available)
        'AnaFlag','StackPixel',...         % FPD segmentations -> pixel combination
        'chi2','chi2Stat',...              % statistics only
        'NonPoissonScaleFactor',1,...      % Non-Poiss. background contribution (1 == no NP background)
        'TwinBias_Q',18573.70,...          % MC endpoint (eV)
        'TwinBias_mnuSq',TwinBias_mnuSq,...             % MC neutrino mass squared (eV^2)
        'DopplerEffectFlag','FSD',...
        'fitter',fitter,...
        'KTFFlag',KTFFlag};
    
    if nStack>1
        % Time bias in (s). 100 times more statistics
        RunAnaArg = {RunAnaArg{:},'TwinBias_Time',7409*nStack};
    end
    
    T = RunAnalysis(RunAnaArg{:});

    if nStack>1
         T.Fit
        FitResult = T.FitResult;
        
        % ScanResults = R.GetAsymFitError('Mode','Uniform',...% equidistant steps
        %                                 'ParScanMax',R.FitResult.err(1)*1.4,...  
        %                                   'SanityPlot','OFF');
        save(savename,'FitResult','R');
        
        if exist('ScanResults','var')
            save(savename,'ScanResults','-append');
        end
    end
end
%% results
if nStack>1
    fprintf('---------------------------------------\n');
    if exist('ScanResults','var')
        fprintf('mNuSq    = %.3f (%.3f +%.3f)eV^2 \n',FitResult.par(1),ScanResults.AsymErr(2),ScanResults.AsymErr(1));
    else
        fprintf('mNuSq    = %.3g (%.3f +%.3f)eV^2 \n',FitResult.par(1),FitResult.errNeg(1),FitResult.errPos(1));
    end
    fprintf('Delta E0 = %.3f (+-%.3f)eV \n',FitResult.par(2)+T.ModelObj.Q_i-18573.70,FitResult.err(2));
    fprintf('N        = %.3f (+-%.4f) \n',FitResult.par(4)+1,FitResult.err(4));
    fprintf('B        = %.1f (+-%.1f) mcps \n',(FitResult.par(3)+T.ModelObj.BKG_RateSec_i)*1e3,1e3*FitResult.err(3));
    fprintf('chi2     = %.1e / %.0f dof \n', FitResult.chi2min,FitResult.dof);
    fprintf('---------------------------------------\n');
    %% Sanity Plot: chi2-curve
    if strcmp(PlotChi2Curve,'ON') && exist('ScanResults','var')
        T.PlotChi2Curve('Parameter','mNu','ScanResult',ScanResults,'FitResult',FitResult);
    end
end
