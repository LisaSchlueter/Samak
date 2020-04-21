%% KATRIN Sterile Neutrino sensitivity - KSN1
%  Nathan Le Guennic - 2020

tic;

%% =====   REAL STAT 40   =====
%  ============================


%% Settings
% Parameters

CL           = 90;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';  % Real
uncertainty  = 'stat';  % syst

d            = 50;                  % Number of dots per decade
eVrange      = 40;                  % eV below the endpoint

savename     = sprintf('coord_%.1deV_%.2s_%.3s_d%.4d.mat',eVrange,datatype,uncertainty,d);


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
        
%% Loop
% Initialisation
%progressbar()
%progressbar(0)
c=0;

diary 'progress.txt'
'STARTING SCAN'
diary off

% Scanning
for m = m4
    % Initializing Newton
    m
    s=1;
    X2_a = fit_chi(m,s,R,eVrange);
    
    if (X2_a>chilvl+err)
        diary 'progress.txt'
        sprintf('Sterile mass : %.1f',m)
        sprintf('\tReducing \\chi^2 ...')
        diary off
        
        % Getting close enough to the value, otherwise minuit is chaotic
        X2_a2=X2_a
        while X2_a2>100 & X2_a2>chilvl
            X2_a=X2_a2;
            s=s/2;
            X2_a2=fit_chi(m,s,R,eVrange);
        end
        
        diary 'progress.txt'
        sprintf('\tStarting Newton convergence\n')
        diary off
        
        % Gradient
        s2   = s-p;
        
        X2_b = fit_chi(m,s2,R,eVrange);

        grad = (X2_a-X2_b)/p;
        beta = min(X2_b,X2_b)-grad*s2 -chilvl;   % Tangeant parameters

        while (X2_a-chilvl>err & X2_b-chilvl>err & s>min_sin2T4 & s<=1)  % Newton loop
            % Gradient computation
            s    = -beta/grad
            X2_a = fit_chi(m,s,R,eVrange)
            
            s2   = s-p
            X2_b = fit_chi(m,s2,R,eVrange)

            % Tangeant parameters
            grad = (X2_a-X2_b)/p;
            beta = min(X2_b,X2_b)-grad*s2 -chilvl;     
        end

        % Actualising plot variables
        if ((abs(X2_a-chilvl)<err|abs(X2_b-chilvl)<err) & s>min_sin2T4 & s<=1)
            sith4_X = [sith4_X,s];
            m4_Y    = [m4_Y,m];
            chi_Z   = [chi_Z,X2_a];
        end
    else
        diary 'progress.txt'
        sprintf('\tOut of bonds')
        diary off
    end
    
    % Progress bar
    c=c+1;
    %progressbar(c/length(m4))
    progress=100*c/length(m4)
    %waitbar(progress,probar,string(progress*100))
    
    diary 'progress.txt'
    Progress=100*c/length(m4)
    diary off
end

t=toc;
diary 'progress.txt'
'SCANNING OVER'
sprintf('Time spent : %.1f hours\n',t/3600)
diary off

%% Datasave
diary 'progress.txt'
'Saving file ...'
diary off

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,savename];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');


diary 'progress.txt'
'==== FINISHED ===='
diary off




%% =====   REAL STAT 90   =====
%  ============================

%% Settings
% Parameters

CL           = 90;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';  % Real
uncertainty  = 'stat';  % syst

d            = 50;                  % Number of dots per decade
eVrange      = 90;                  % eV below the endpoint

savename     = sprintf('coord_%.1deV_%.2s_%.3s_d%.4d.mat',eVrange,datatype,uncertainty,d);


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
        
%% Loop
% Initialisation
%progressbar()
%progressbar(0)
c=0;

diary 'progress.txt'
'STARTING SCAN'
diary off

% Scanning
for m = m4
    % Initializing Newton
    m
    s=1;
    X2_a = fit_chi(m,s,R,eVrange);
    
    if (X2_a>chilvl+err)
        diary 'progress.txt'
        sprintf('Sterile mass : %.1f',m)
        sprintf('\tReducing \\chi^2 ...')
        diary off
        
        % Getting close enough to the value, otherwise minuit is chaotic
        X2_a2=X2_a
        while X2_a2>100 & X2_a2>chilvl
            X2_a=X2_a2;
            s=s/2;
            X2_a2=fit_chi(m,s,R,eVrange);
        end
        
        diary 'progress.txt'
        sprintf('\tStarting Newton convergence\n')
        diary off
        
        % Gradient
        s2   = s-p;
        
        X2_b = fit_chi(m,s2,R,eVrange);

        grad = (X2_a-X2_b)/p;
        beta = min(X2_b,X2_b)-grad*s2 -chilvl;   % Tangeant parameters

        while (X2_a-chilvl>err & X2_b-chilvl>err & s>min_sin2T4 & s<=1)  % Newton loop
            % Gradient computation
            s    = -beta/grad
            X2_a = fit_chi(m,s,R,eVrange)
            
            s2   = s-p
            X2_b = fit_chi(m,s2,R,eVrange)

            % Tangeant parameters
            grad = (X2_a-X2_b)/p;
            beta = min(X2_b,X2_b)-grad*s2 -chilvl;     
        end

        % Actualising plot variables
        if ((abs(X2_a-chilvl)<err|abs(X2_b-chilvl)<err) & s>min_sin2T4 & s<=1)
            sith4_X = [sith4_X,s];
            m4_Y    = [m4_Y,m];
            chi_Z   = [chi_Z,X2_a];
        end
    else
        diary 'progress.txt'
        sprintf('\tOut of bonds')
        diary off
    end
    
    % Progress bar
    c=c+1;
    %progressbar(c/length(m4))
    progress=100*c/length(m4)
    %waitbar(progress,probar,string(progress*100))
    
    diary 'progress.txt'
    Progress=100*c/length(m4)
    diary off
end

t=toc;
diary 'progress.txt'
'SCANNING OVER'
sprintf('Time spent : %.1f hours\n',t/3600)
diary off

%% Datasave
diary 'progress.txt'
'Saving file ...'
diary off

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,savename];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');


diary 'progress.txt'
'==== FINISHED ===='
diary off




%% =====   REAL SYST 40   =====
%  ============================

%% Settings
% Parameters

CL           = 90;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';  % Real
uncertainty  = 'syst';  % syst

d            = 50;                  % Number of dots per decade
eVrange      = 40;                  % eV below the endpoint

savename     = sprintf('coord_%.1deV_%.2s_%.3s_d%.4d.mat',eVrange,datatype,uncertainty,d);


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
        
%% Loop
% Initialisation
%progressbar()
%progressbar(0)
c=0;

diary 'progress.txt'
'STARTING SCAN'
diary off

% Scanning
for m = m4
    % Initializing Newton
    m
    s=1;
    X2_a = fit_chi(m,s,R,eVrange);
    
    if (X2_a>chilvl+err)
        diary 'progress.txt'
        sprintf('Sterile mass : %.1f',m)
        sprintf('\tReducing \\chi^2 ...')
        diary off
        
        % Getting close enough to the value, otherwise minuit is chaotic
        X2_a2=X2_a
        while X2_a2>100 & X2_a2>chilvl
            X2_a=X2_a2;
            s=s/2;
            X2_a2=fit_chi(m,s,R,eVrange);
        end
        
        diary 'progress.txt'
        sprintf('\tStarting Newton convergence\n')
        diary off
        
        % Gradient
        s2   = s-p;
        
        X2_b = fit_chi(m,s2,R,eVrange);

        grad = (X2_a-X2_b)/p;
        beta = min(X2_b,X2_b)-grad*s2 -chilvl;   % Tangeant parameters

        while (X2_a-chilvl>err & X2_b-chilvl>err & s>min_sin2T4 & s<=1)  % Newton loop
            % Gradient computation
            s    = -beta/grad
            X2_a = fit_chi(m,s,R,eVrange)
            
            s2   = s-p
            X2_b = fit_chi(m,s2,R,eVrange)

            % Tangeant parameters
            grad = (X2_a-X2_b)/p;
            beta = min(X2_b,X2_b)-grad*s2 -chilvl;     
        end

        % Actualising plot variables
        if ((abs(X2_a-chilvl)<err|abs(X2_b-chilvl)<err) & s>min_sin2T4 & s<=1)
            sith4_X = [sith4_X,s];
            m4_Y    = [m4_Y,m];
            chi_Z   = [chi_Z,X2_a];
        end
    else
        diary 'progress.txt'
        sprintf('\tOut of bonds')
        diary off
    end
    
    % Progress bar
    c=c+1;
    %progressbar(c/length(m4))
    progress=100*c/length(m4)
    %waitbar(progress,probar,string(progress*100))
    
    diary 'progress.txt'
    Progress=100*c/length(m4)
    diary off
end

t=toc;
diary 'progress.txt'
'SCANNING OVER'
sprintf('Time spent : %.1f hours\n',t/3600)
diary off

%% Datasave
diary 'progress.txt'
'Saving file ...'
diary off

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,savename];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');


diary 'progress.txt'
'==== FINISHED ===='
diary off


%% =====   REAL SYST 90 90%  =====
%  ==============================

%% Settings
% Parameters

CL           = 90;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';  % Real
uncertainty  = 'syst';  % syst

d            = 50;                  % Number of dots per decade
eVrange      = 90;                  % eV below the endpoint

savename     = sprintf('coord_%.1deV_%.2s_%.3s_d%.4d.mat',eVrange,datatype,uncertainty,d);


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
        
%% Loop
% Initialisation
%progressbar()
%progressbar(0)
c=0;

diary 'progress.txt'
'STARTING SCAN'
diary off

% Scanning
for m = m4
    % Initializing Newton
    m
    s=1;
    X2_a = fit_chi(m,s,R,eVrange);
    
    if (X2_a>chilvl+err)
        diary 'progress.txt'
        sprintf('Sterile mass : %.1f',m)
        sprintf('\tReducing \\chi^2 ...')
        diary off
        
        % Getting close enough to the value, otherwise minuit is chaotic
        X2_a2=X2_a
        while X2_a2>100 & X2_a2>chilvl
            X2_a=X2_a2;
            s=s/2;
            X2_a2=fit_chi(m,s,R,eVrange);
        end
        
        diary 'progress.txt'
        sprintf('\tStarting Newton convergence\n')
        diary off
        
        % Gradient
        s2   = s-p;
        
        X2_b = fit_chi(m,s2,R,eVrange);

        grad = (X2_a-X2_b)/p;
        beta = min(X2_b,X2_b)-grad*s2 -chilvl;   % Tangeant parameters

        while (X2_a-chilvl>err & X2_b-chilvl>err & s>min_sin2T4 & s<=1)  % Newton loop
            % Gradient computation
            s    = -beta/grad
            X2_a = fit_chi(m,s,R,eVrange)
            
            s2   = s-p
            X2_b = fit_chi(m,s2,R,eVrange)

            % Tangeant parameters
            grad = (X2_a-X2_b)/p;
            beta = min(X2_b,X2_b)-grad*s2 -chilvl;     
        end

        % Actualising plot variables
        if ((abs(X2_a-chilvl)<err|abs(X2_b-chilvl)<err) & s>min_sin2T4 & s<=1)
            sith4_X = [sith4_X,s];
            m4_Y    = [m4_Y,m];
            chi_Z   = [chi_Z,X2_a];
        end
    else
        diary 'progress.txt'
        sprintf('\tOut of bonds')
        diary off
    end
    
    % Progress bar
    c=c+1;
    %progressbar(c/length(m4))
    progress=100*c/length(m4)
    %waitbar(progress,probar,string(progress*100))
    
    diary 'progress.txt'
    Progress=100*c/length(m4)
    diary off
end

t=toc;
diary 'progress.txt'
'SCANNING OVER'
sprintf('Time spent : %.1f hours\n',t/3600)
diary off

%% Datasave
diary 'progress.txt'
'Saving file ...'
diary off

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,savename];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');


diary 'progress.txt'
'==== FINISHED ===='
diary off



%% =====   REAL SYST 90 95%  =====
%  ==============================

%% Settings
% Parameters

CL           = 95;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';  % Real
uncertainty  = 'syst';  % syst

d            = 50;                  % Number of dots per decade
eVrange      = 90;                  % eV below the endpoint

savename     = sprintf('coord_%.1deV_%.2s_%.3s_d%.4d_95%.mat',eVrange,datatype,uncertainty,d);


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
        
%% Loop
% Initialisation
%progressbar()
%progressbar(0)
c=0;

diary 'progress.txt'
'STARTING SCAN'
diary off

% Scanning
for m = m4
    % Initializing Newton
    m
    s=1;
    X2_a = fit_chi(m,s,R,eVrange);
    
    if (X2_a>chilvl+err)
        diary 'progress.txt'
        sprintf('Sterile mass : %.1f',m)
        sprintf('\tReducing \\chi^2 ...')
        diary off
        
        % Getting close enough to the value, otherwise minuit is chaotic
        X2_a2=X2_a
        while X2_a2>100 & X2_a2>chilvl
            X2_a=X2_a2;
            s=s/2;
            X2_a2=fit_chi(m,s,R,eVrange);
        end
        
        diary 'progress.txt'
        sprintf('\tStarting Newton convergence\n')
        diary off
        
        % Gradient
        s2   = s-p;
        
        X2_b = fit_chi(m,s2,R,eVrange);

        grad = (X2_a-X2_b)/p;
        beta = min(X2_b,X2_b)-grad*s2 -chilvl;   % Tangeant parameters

        while (X2_a-chilvl>err & X2_b-chilvl>err & s>min_sin2T4 & s<=1)  % Newton loop
            % Gradient computation
            s    = -beta/grad
            X2_a = fit_chi(m,s,R,eVrange)
            
            s2   = s-p
            X2_b = fit_chi(m,s2,R,eVrange)

            % Tangeant parameters
            grad = (X2_a-X2_b)/p;
            beta = min(X2_b,X2_b)-grad*s2 -chilvl;     
        end

        % Actualising plot variables
        if ((abs(X2_a-chilvl)<err|abs(X2_b-chilvl)<err) & s>min_sin2T4 & s<=1)
            sith4_X = [sith4_X,s];
            m4_Y    = [m4_Y,m];
            chi_Z   = [chi_Z,X2_a];
        end
    else
        diary 'progress.txt'
        sprintf('\tOut of bonds')
        diary off
    end
    
    % Progress bar
    c=c+1;
    %progressbar(c/length(m4))
    progress=100*c/length(m4)
    %waitbar(progress,probar,string(progress*100))
    
    diary 'progress.txt'
    Progress=100*c/length(m4)
    diary off
end

t=toc;
diary 'progress.txt'
'SCANNING OVER'
sprintf('Time spent : %.1f hours\n',t/3600)
diary off

%% Datasave
diary 'progress.txt'
'Saving file ...'
diary off

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,savename];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');


diary 'progress.txt'
'==== FINISHED ===='
diary off




%% =====   REAL SYST 90 99%  =====
%  ==============================

%% Settings
% Parameters

CL           = 99;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';  % Real
uncertainty  = 'syst';  % syst

d            = 50;                  % Number of dots per decade
eVrange      = 90;                  % eV below the endpoint

savename     = sprintf('coord_%.1deV_%.2s_%.3s_d%.4d_99%.mat',eVrange,datatype,uncertainty,d);


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
        
%% Loop
% Initialisation
%progressbar()
%progressbar(0)
c=0;

diary 'progress.txt'
'STARTING SCAN'
diary off

% Scanning
for m = m4
    % Initializing Newton
    m
    s=1;
    X2_a = fit_chi(m,s,R,eVrange);
    
    if (X2_a>chilvl+err)
        diary 'progress.txt'
        sprintf('Sterile mass : %.1f',m)
        sprintf('\tReducing \\chi^2 ...')
        diary off
        
        % Getting close enough to the value, otherwise minuit is chaotic
        X2_a2=X2_a
        while X2_a2>100 & X2_a2>chilvl
            X2_a=X2_a2;
            s=s/2;
            X2_a2=fit_chi(m,s,R,eVrange);
        end
        
        diary 'progress.txt'
        sprintf('\tStarting Newton convergence\n')
        diary off
        
        % Gradient
        s2   = s-p;
        
        X2_b = fit_chi(m,s2,R,eVrange);

        grad = (X2_a-X2_b)/p;
        beta = min(X2_b,X2_b)-grad*s2 -chilvl;   % Tangeant parameters

        while (X2_a-chilvl>err & X2_b-chilvl>err & s>min_sin2T4 & s<=1)  % Newton loop
            % Gradient computation
            s    = -beta/grad
            X2_a = fit_chi(m,s,R,eVrange)
            
            s2   = s-p
            X2_b = fit_chi(m,s2,R,eVrange)

            % Tangeant parameters
            grad = (X2_a-X2_b)/p;
            beta = min(X2_b,X2_b)-grad*s2 -chilvl;     
        end

        % Actualising plot variables
        if ((abs(X2_a-chilvl)<err|abs(X2_b-chilvl)<err) & s>min_sin2T4 & s<=1)
            sith4_X = [sith4_X,s];
            m4_Y    = [m4_Y,m];
            chi_Z   = [chi_Z,X2_a];
        end
    else
        diary 'progress.txt'
        sprintf('\tOut of bonds')
        diary off
    end
    
    % Progress bar
    c=c+1;
    %progressbar(c/length(m4))
    progress=100*c/length(m4)
    %waitbar(progress,probar,string(progress*100))
    
    diary 'progress.txt'
    Progress=100*c/length(m4)
    diary off
end

t=toc;
diary 'progress.txt'
'SCANNING OVER'
sprintf('Time spent : %.1f hours\n',t/3600)
diary off

%% Datasave
diary 'progress.txt'
'Saving file ...'
diary off

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,savename];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');


diary 'progress.txt'
'==== FINISHED ===='
diary off


%%
'================================'
'===    EVERYTHING FINISHED   ==='
'================================'