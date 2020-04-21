%%
% A tool to scan only one choosen line

%% Sterile mass
% m = 686.649;
m = 68.6649;

%% Parameters
CL           = 90;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';  % Twin
uncertainty  = 'stat';  % syst

eVrange      = 40;                  % eV below the endpoint

min_sin2T4   = 0.00001;              % Lower bound for sin2(th4) for the plot

p   = 0.00001;                       % Newton gradient step
err = 0.01;                         % Newton convergence tolerance

% Confidence level
switch CL
    case 90
        chilvl = 4.61;
    case 95
        chilvl = 5.99;
    case 99
        chilvl = 9.21;
end

% Uncertainties
switch uncertainty
    case 'stat'
        chi2_type = 'chi2Stat';
    case 'syst'
        chi2_type = 'chi2CMShape';
end

%% Constructor Initialisation
switch datatype
    case 'Real'
        sibille = 'SibilleFull';
        chilvl  = chilvl+22.7;      % Delta chi2
    case 'Twin'
        sibille = 'Sibille0p5eV';
end
R = MultiRunAnalysis('RunList','KNM1',...
            'chi2',chi2_type,...
            'DataType',datatype,...
            'fixPar','E0 Norm Bkg',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag',sibille,...
            'ELossFlag','KatrinT2',...
            'SysBudget',22);
        
% Initializing Newton
s=1;
X2_a = fit_chi(m,s,R,eVrange);

fprintf('\nSterile mass : %.3f',m)

if (X2_a>chilvl+err)
    fprintf('\n\tReducing \\chi^2 ...\n')

    % Getting close enough to the value, otherwise minuit is chaotic
    X2_a2=X2_a
    while X2_a2>200
        s=s/1.3;
        X2_a2=fit_chi(m,s,R,eVrange);
    end

    fprintf('\tStarting Newton convergence\n')

    % Gradient
    X2_a = X2_a2
    s2   = s-p;
    X2_b = fit_chi(m,s2,R,eVrange)

    if abs(X2_a-X2_b)/X2_b<0.1          % Prevent Chaos

        grad = (X2_a-X2_b)/p;
        beta = X2_a-grad*s -chilvl;     % Tangeant parameters

        while (X2_a-chilvl>err & X2_b-chilvl>err & s>min_sin2T4 & s<=1)  % Newton loop
            % Gradient computation
            s    = -beta/grad;
            X2_a = fit_chi(m,s,R,eVrange)
            
            if X2_a>5*X2_b              % Security
                X2_a=0;
            end
            
            s2   = s-p;
            X2_b = fit_chi(m,s2,R,eVrange)

            % Tangeant parameters
            grad = (X2_a-X2_b)/p;
            beta = X2_a-grad*s -chilvl;     
        end
    end
else
    fprintf('\n\tOut of bonds\n')
end

disp('OVER');