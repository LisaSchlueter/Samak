%% Exmaple of a linear fit
% This example was designed to follow a similar design schema as given by more complex
% minimization problems for analyzing the tritium beta decay model. 
% Here a simple linear graph class is used: 'LinearSim.m'. The simulation class defines a method 
% 'compute(vargargin)' which returns an array of values and accepts two parameters in the form of 
% ''name',value' pairs: 'slope' and 'offset'. 
% The minimization function comparing the theory and data is defined in 'tools/fitter/minimize.m'. 
% This function is minimized using the minimization program MINUIT ('tools/fitter/fminuit.m').
%
% 2017 Nov - Marc Korzeczek

%% Prepare script environment
clearvars; delete(findall(0,'Type','figure'));

addpath('../../tools/fitter');  %for fminuit() and minimize.chi2()

%% Prepare data and theory
sim     = LinearSim(); %value-class constructor

data    = sim.compute('slope',10,'offset',500,'poisson',true); %return array of data
sigma   = sqrt(data);                                          %array of uncertainties
theoFun = @sim.compute;                                        %function handle for array of theory
fitPar  = {'slope','offset'};                                  %parameters fitted in 'theoFun'
auxChi2Fun=@(slope,offset) (slope-10)^2/sqrt(10)+(offset-500)/sqrt(500);

% Define the minimizer function passed to minuit and some useful arguments
funName   = 'minimize.chi2';             % name of function that will be minimized
funInitGuess = [.3, 0];                  % initial guess of fit parameters
funArg  = {data,sigma,theoFun,fitPar};   % cell of mandatory arguments for 'FunName'
funArg  = [funArg,{auxChi2Fun}]; % if needed add auxiliary chi2-terms 

% Call minimization program minuit
% - Minuit will internally minimize 'funName(funInitGuess, funArg)'
% - Check Minuit documentation for more information about arguments
minuitNames = cell2mat(strcat(fitPar,{' '}));
minuitArg = {'-n',minuitNames,'-c','set pri 0; min'};

out = fminuit( funName, funInitGuess, funArg, minuitArg{:});


%% Plot original data and fit-result
figure(1);
x=.1:.1:10;
plot(x,data,'b+','MarkerSize',2); hold on;

tmp=num2cell(out); fitArgIn={fitPar{:}; tmp{:}};
plot(x,theoFun(fitArgIn{:})); hold off;
