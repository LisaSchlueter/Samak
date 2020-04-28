% test parallel grid search for sterile analysis
% Lisa, April 2020
savedir = [getenv('SamakPath'),'ksn1ana/LisaSterile/results/'];
MakeDir(savedir);
range = 40;
nGridSteps = 50;
NonPoissonScaleFactor=1;
chi2 = 'chi2Stat';
DataType = 'Real';
freePar = 'E0 Bkg Norm';
RunList = 'KNM1';
savefile = sprintf('%sSterileTestPar_%s_%s_%s_%.0feVrange_%.0fGridSteps.mat',...
    savedir,RunList,DataType,strrep(freePar,' ',''),range,nGridSteps);
if exist(savefile,'file')
else
    RunAnaArg = {'RunList',RunList,...
        'fixPar',freePar,...
        'DataType',DataType,...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'AnaFlag','StackPixel',...
        'chi2',chi2,...
        'ROIFlag','Default',...
        'SynchrotronFlag','ON',...
        'AngularTFFlag','OFF',...
        'ISCSFlag','Edep',...
        'NonPoissonScaleFactor',NonPoissonScaleFactor};
    T = MultiRunAnalysis(RunAnaArg{:});
    T.exclDataStart = T.GetexclDataStart(range);
    
    T.Fit;
    FitResults_ref = T.FitResult;
    chi2_ref       = T.FitResult.chi2min;
    %% define grid
    mnu4Sq = repmat(linspace(0,50^2,nGridSteps)',1,nGridSteps);
    sin2T4 = repmat(linspace(0.001,1,nGridSteps),nGridSteps,1); %logspace(-3,0,nGridSteps)
    
    mnu4Sq_Grid = reshape(mnu4Sq',nGridSteps.*nGridSteps,1);
    sin2T4_Grid = reshape(sin2T4',nGridSteps.*nGridSteps,1); 
    % make copy of models for parallel computing
    D = copy(repmat(T,nGridSteps.*nGridSteps,1));
   % D = gpuArray(D);
    %arrayfun(@(x) x.SimulateStackRuns,D'); % get model object
   % arrayfun(@(x,mu,sin) x.ModelObj.SetFitBiasSterile(mu,sin),D',mnu4Sq,sin2T4); % set sterile parameters
   
  %  arrayfun(@(x,mu,sin) x.ModelObj.SetFitBiasSterile(mu,sin),D,mnu4Sq_Grid,sin2T4_Grid); % set sterile parameters
   
   %%
   D             = reshape(D,numel(D),1);
   chi2Grid       = zeros(nGridSteps*nGridSteps,1);
   FitResultsGrid = cell(nGridSteps*nGridSteps,1);
   
   parfor i= 1:(nGridSteps*nGridSteps)
       D(i).SimulateStackRuns;
       D(i).ModelObj.SetFitBiasSterile(mnu4Sq_Grid(i),sin2T4_Grid(i));
       D(i).Fit
       chi2Grid(i) = D(i).FitResult.chi2min;
       FitResultsGrid{i} = D(i).FitResult;
   end
  
  D              = reshape(D,nGridSteps,nGridSteps);
  chi2Grid       = reshape(chi2Grid,nGridSteps,nGridSteps);
  FitResultsGrid = reshape(FitResultsGrid,nGridSteps,nGridSteps);
    %%
%     chi2Grid        = arrayfun(@(x) x.FitResult.chi2min,D);
%     FitResultsGrid  = arrayfun(@(x) x.FitResult,D);
  % save(savefile,'chi2_ref','FitResults_ref','RunAnaArg',...
  %     'chi2Grid','mnu4Sq','sin2T4','FitResultsGrid');
  % fprintf('save file to %s \n',savefile)
end

% %% get contour at X sigma
% CL = 90;
% % Confidence level
% switch CL
%     case 90
%         DeltaChi2 = 4.61;
%     case 95
%         DeltaChi2 = 5.99;
%     case 99
%         DeltaChi2 = 9.21;
% end
% 
% sin2T4_contour = zeros(nGridSteps,1);
% mnu4Sq_contour = mnu4Sq(:,1);
% 
% for i=1:numel(nGridSteps)
%     sin2T4_contour(i) = interp1(chi2Grid(:,i)-chi2_ref,mnu4Sq_contour,DeltaChi2,'spline');
% end
