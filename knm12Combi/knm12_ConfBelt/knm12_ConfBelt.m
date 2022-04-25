% confidence belt for combined analysis via summation of chi2 profiles
% method:
% 1. load standalone profiles (asimov)
% 2. combine

savedir = [getenv('SamakPath'),'knm12Combi/knm12_ConfBelt/results/'];
savename = [savedir,'knm12_ConfBelt.mat'];

Sensitivity = 'ON'; % OFF= show best fit, ON = show sensitivity only
SavePlot    = 'OFF';
range = 40;
mNuSq = sort([0:0.1:0.6,0.8,0.95,0.75,0.65,0.15,1.35]);


if exist(savename,'file')
    load(savename);
else
    
    %%  set up knm2 model
    mNuSq2_bf = 0.279; % uniform
    RunAnaArg_KNM2 = {'RunList','KNM2_Prompt',...
        'fixPar','mNu E0 Norm Bkg',...%free par
        'SysBudget',40,...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'RadiativeFlag','ON',...
        'FSDFlag','KNM2',...
        'ELossFlag','KatrinT2A20',...
        'AnaFlag','StackPixel',...
        'chi2','chi2CMShape',...
        'NonPoissonScaleFactor',1.112,...
        'FSD_Sigma',sqrt(0.0124+0.0025),...
        'TwinBias_FSDSigma',sqrt(0.0124+0.0025),...
        'TwinBias_Q',18573.7,...
        'PullFlag',99,...;%99 = no pull
        'BKG_PtSlope', 3*1e-06,...
        'TwinBias_BKG_PtSlope', 3*1e-06,...
        'DopplerEffectFlag','FSD'};
    
    M2 = MultiRunAnalysis(RunAnaArg_KNM2{:});
    M2.exclDataStart = M2.GetexclDataStart(range);
    S2 = RunSensitivity('RunAnaObj',M2);
    S2.ConfLevel = 0.9; % confidence level (0==1 sigma)
    
    %% knm1
    RunAnaArg_KNM1 = {'RunList','KNM1',...
        'fixPar','mNu E0 Norm Bkg',...%free par
        'exclDataStart',13,...
        'SysBudget',22,...
        'FSDFlag','Sibille0p5eV',...
        'fitter','minuit',...
        'minuitOpt','min;migrad',...
        'NonPoissonScaleFactor',1.064,...
        'RadiativeFlag','ON',...
        'ELossFlag','KatrinT2',...
        'chi2','chi2CMShape',...
        'AngularTF','OFF'};
    % set up model
    M1 = MultiRunAnalysis(RunAnaArg_KNM1{:});
    S1 = RunSensitivity('RunAnaObj',M1);
    
    %% Compute and Plot ConfidencE Belt
    S1.ConfLevel = 0.9; % confidence level (0==1 sigma)
    save(savename,'S1','S2');
end
%% combine
Mode        = 'FC';  % FC = Feldman Cousin, LT = Lokov Tkachov
mNuSq = sort([0:0.1:0.6,0.95,0.8,0.75,0.65,0.15,1.35]);

for i=1:numel(mNuSq)
    % Get Delta Chi2 Curve
    [mNuSq_x1,DeltaChi21,Chi2True1,Chi2Best1] = S1.FC_ComputeDeltaChi2LookupTables('mNuSq_t',mNuSq(i),'nSamples',300);
    [mNuSq_x2,DeltaChi22,Chi2True2,Chi2Best2] = S2.FC_ComputeDeltaChi2LookupTables('mNuSq_t',mNuSq(i),'nSamples',300);
    
    % interpolate
    xmin  = max([min(mNuSq_x1) min(mNuSq_x1)]);
    xmax  = min([max(mNuSq_x1) max(mNuSq_x1)]);
    mNuSq_x   = linspace(xmin,xmax,300);
    DeltaChi2 = interp1(mNuSq_x1,DeltaChi21,mNuSq_x,'lin')+interp1(mNuSq_x2,DeltaChi22,mNuSq_x,'lin');
    Chi2True = Chi2True1+Chi2True2;
    Prob_tmp = exp(-Chi2True./2);
    Prob = Prob_tmp./simpsons(mNuSq_x,Prob_tmp); % normalization
    CumProb = GetCDF(mNuSq_x,Prob);%CumProb = cumsum(Prob); CumProb = CumProb./max(CumProb);
    
    switch Mode
        case 'FC'
            % solve (numerically) after x1 and x2
            options = optimoptions(@fsolve,'Display','off');
            interpStyle = 'lin';
            if mNuSq(i)~=0
                a =  @(x) (interp1(mNuSq_x,CumProb,x(2),interpStyle)-interp1(mNuSq_x,CumProb,x(1),interpStyle))-S2.ConfLevel;
                b =  @(x)  (interp1(mNuSq_x,DeltaChi2,x(1),interpStyle)-interp1(mNuSq_x,DeltaChi2,x(2),interpStyle));
                fun = @(x) [a(x),b(x)];
                x0 = [-1,1];
                x1x2 = fsolve(fun,x0,options);
                S2.FC_x1(i) = x1x2(1);
                S2.FC_x2(i) = x1x2(2);
                S2.FC_DeltaChi2C(i) = interp1(mNuSq_x,DeltaChi2,S2.FC_x1(i),interpStyle);
            elseif mNuSq(i)==0 % Delta Chi2 is always zero for negative masses
                fun =  @(x) interp1(mNuSq_x,CumProb,x,interpStyle)-S2.ConfLevel;
                x1x2 = fsolve(fun,1,options);
                S2.FC_x1(i) = NaN;
                S2.FC_x2(i) = x1x2;
                S2.FC_DeltaChi2C(i) = interp1(mNuSq_x,DeltaChi2,S2.FC_x2(i),interpStyle);
            end
        case 'LT'
            % find acceptance region
            [CumProb,ia,~] = unique(CumProb);
            S2.FC_x1(i) = interp1(CumProb,mNuSq_x(ia),(1-S2.ConfLevel)/2,'spline');
            
            if i~=1
                if S2.FC_x1(i)-S2.FC_x1(i-1)<0
                    % should always be more positive! try piecewise
                    % cubic interpolation instead
                    S2.FC_x1(i) = interp1(CumProb,mNuSq_x(ia),(1-S2.ConfLevel)/2,'pchip');
                end
            end
            
            if S2.FC_x1(i)>0  % --> twosided
                S2.FC_x2(i)= interp1(CumProb,mNuSq_x(ia),(1+S2.ConfLevel)/2,'spline');
            elseif S2.FC_x1(i)<0 % --> one-sided
                S2.FC_x2(i) = interp1(CumProb,mNuSq_x(ia),S2.ConfLevel,'spline');
            end
    end
end

S2.FC_mNuSqTrue = mNuSq';

%%
switch Mode
    case 'FC'
        S2.PlotFCBelt('HoldOn','OFF','Sensitivity','OFF',...
            'SavePlot','OFF','XLim',[-1,1],...
            'Style','Pretty','mNuSq_bf',0.10);
    case 'LT'
        S2.PlotFCBelt('Lokov','ON','Sensitivity',Sensitivity,'SavePlot',SavePlot,...
            'Style','Pretty','XLim',[-1.25,1.25],'mNuSq_bf',0.10);
end
