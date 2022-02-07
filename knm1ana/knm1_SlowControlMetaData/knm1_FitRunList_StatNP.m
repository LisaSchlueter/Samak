% fit knm1 run list (again) with stat + NonPoisson

% KNM1 Final fit Results
% Uniform Fit
% Golden Run List
% Golden Pixel List

savedir = [getenv('SamakPath'),'knm1ana/knm1_SlowControlMetaData/results/'];
savefile = sprintf('%sknm1_RunwiseFits_NP.mat',savedir);

if exist(savefile,'file')
    load(savefile);
else
    
    %% settings
    RunList = 'KNM1';
    range   = 40;         % 40eV range = 27 subruns
    % Init Model Object and covariance matrix object
    R = MultiRunAnalysis('RunList',RunList,...
        'chi2','chi2Stat',...
        'DataType','Real',...
        'fixPar','E0 Norm Bkg',... free parameter
        'RadiativeFlag','ON',...
        'NonPoissonScaleFactor',1.064,...
        'minuitOpt','min ; minos',...
        'FSDFlag','Sibille0p5eV',...
        'ELossFlag','KatrinT2',...
        'AngularTFFlag','OFF');
    R.exclDataStart = R.GetexclDataStart(range);
    FitResults = R.FitRunList;
    save(savefile,'R','FitResults');
end


