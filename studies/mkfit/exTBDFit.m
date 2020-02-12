%% Example of TBD Fit
%
%  Fit For Neutrino Mass
%  Fit Class By. Marc Korzeczek
%  Last Update: 27/02/2018
%

%% Prepare script environment
clearvars; delete(findall(0,'Type','figure'));

addpath('../../simulation/tritium');%for exLinearSim()
addpath('../../fitting');   %for fminuit() and minimize.chi2()
addpath('../../tools');     %for figureFormat()

%% Prepare data and theory
%model
sim  = TBD('TimeSec',86400*1000);         
sim.ComputeTBDDS; sim.ComputeTBDIS; 

%data
D    = TBD('TimeSec',86400*1000,'mnuSq_i',sqrt(0.42));         
D.ComputeTBDDS; D.ComputeTBDIS; D.AddStatFluctTBDIS;
data    = D.TBDIS;              %array of data
sigma   = sim.TBDISE;             %array of uncertainties

% Define the minimizer function passed to minuit and some useful arguments
theoFun = @sim.FitCalculateTBDIS;%function handle for array of theory
fitPar  = {'mSq_bias','E0_bias','B_bias','N_bias'};%parameters fitted in 'theoFun'
funName   = 'minimize.chi2Gauss';                  % name of function that will be minimized
funInitGuess = [ 0, 0, 0, 0];                  % initial guess of fit parameters

% Call minimization program minuit
% - Minuit will internally minimize 'funName(funInitGuess, funArg)'
minuitNames = cell2mat(strcat(fitPar,{' '}));
minuitArg = {'-n',minuitNames,'-c','set pri 1; minos'};

nfit = 5;
outstat = zeros(nfit,numel(funInitGuess)+1);
for i=1:1:nfit
timerVal = cputime;
% data 
D.SetFitBias(0); D.ComputeTBDDS; D.ComputeTBDIS; D.AddStatFluctTBDIS;
data    = D.TBDIS;              %array of data
sigma   = D.TBDISE;             %array of uncertainties

% Minimization
funArg  = {data,sigma,theoFun,fitPar};% cell of mandatory arguments for 'FunName'
[out, err, chi2min, errmat] = fminuit( funName, funInitGuess, funArg, minuitArg{:});
timerVal = cputime - timerVal;

%% Plot original data and fit-result
figure(1);
x=sim.qU;
tmp=num2cell(out); fitArgIn={fitPar{:}; tmp{:}};
plot=errorbar(x,(data-theoFun(fitArgIn{:})),sqrt(data),'b+','MarkerSize',2);
strleg = sprintf('Fit %g (%g/%g/%g) - \\chi^2=%g/%g- CPU=%g',i,out(1),out(2),out(3),chi2min,numel(x)-numel(out),timerVal);
legend(plot,strleg);
%plot(x,theoFun(fitArgIn{:})); 
outstat(i,:) = [out chi2min];
end

figure(2)
subplot(2,2,1)
nhist(outstat(:,1));
subplot(2,2,2)
nhist(outstat(:,2));
subplot(2,2,3)
nhist(outstat(:,3));
subplot(2,2,4)
nhist(outstat(:,4));
