%% KATRIN Sterile Neutrino sensitivity - KSN1
%  Sanity checks
%  Nathan Le Guennic - 2020

tic;
%% Settings
% Parameters
savename     = 'coord_0to4_d400_90eV_syst.mat';
CL           = 90;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';
uncertainty  = 'chi2CMShape'; % chi2Stat

d            = 200;                 % Number of dots per decade
eVrange      = 90;                  % eV below the endpoint


start_decade = 0;
stop_decade  = 4;
min_sin2T4   = 0.0001;              % Lower bound for sin2(th4) for the plot


p   = 0.0000001;                    % Newton gradient step
err = 0.01;                         % Newton convergence tolerance


% Variables
m4_Y    = [];
sith4_X = [];                       % Plot variables
chi_Z   = [];

% Scanning domain
m4      = [];
k=1;
for n=start_decade:stop_decade-2
    m4=[m4,logspace(n,n+1,d)];
    m4(length(m4))=[];
    k=k+1;
end % Deleting the middle points
m4 = [m4,logspace(stop_decade-1,stop_decade,d)];


% Confidence level
switch CL
    case 90
        chilvl = 4.61;
    case 95
        chilvl = 5.99;
    case 99
        chilvl = 9.21;
end

%% Constructor Initialisation
switch datatype
    case 'Real'
        sibille = 'SibilleFull';
    case 'Twin'
        sibille = 'Sibille0p5eV';
end

R = MultiRunAnalysis('RunList','KNM1',...
            'chi2',uncertainty,...
            'DataType',datatype,...
            'fixPar','E0 Norm Bkg',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',1.064,...
            'minuitOpt','min ; migrad',...
            'FSDFlag',sibille,...
            'ELossFlag','KatrinT2',...
            'SysBudget',22);
        
%% Loop
% Initialisation
%progressbar()
%progressbar(0)
c=0;

% Scanning
for m = m4
    % Initializing Newton
    m
    s=1;
    X2_a = fit_chi(m,s,R,eVrange)
    
    if (X2_a>chilvl+err)
        % Gradient
        s = s-p;
        X2_b = fit_chi(m,s,R,eVrange)
        
        if abs(X2_b-X2_a)/min(X2_a,X2_b)<1        % Preventing some computational errors
            
            grad = (X2_a-X2_b)/p;
            beta = X2_b-grad*s -chilvl;   % Tangeant parameters

            while (X2_a-chilvl>err & X2_b-chilvl>err & s>min_sin2T4 & s<=1)  % Newton loop
                % Gradient computation
                s = -beta/grad;
                X2_a = fit_chi(m,s,R,eVrange)

                s = s-p;
                X2_b = fit_chi(m,s,R,eVrange)

                % Tangeant parameters
                grad = (X2_a-X2_b)/p;
                beta = X2_b-grad*s -chilvl;     
            end

            % Actualising plot variables
            if ((abs(X2_a-chilvl)<err|abs(X2_b-chilvl)<err) & s>min_sin2T4 & s<=1)
                sith4_X = [sith4_X,s];
                m4_Y    = [m4_Y,m];
                chi_Z   = [chi_Z,X2_a];
            end
        end
    end
    
    % Progress
    c=c+1;
    progress=100*c/length(m4)
    diary 'progress.txt'
    progress=100*c/length(m4)
    diary off
end

t=toc;
'SCANNING OVER'
t/3600

%% Datasave
filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,savename];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');

'FINISHED'
