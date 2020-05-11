%% KATRIN Sterile Neutrino sensitivity - KSN1
%  Nathan Le Guennic - 2020

tic;
%% Settings
% Parameters

CL           = 95;                  % Confidence level 90% - 95% - 99%
datatype     = 'Real';  % Real
uncertainty  = 'syst';  % syst
NPfactor     = 1.064;
d            = 10;                  % Number of dots per decade
eVrange      = 40;                  % eV below the endpoint

% Name for the datafile
%savename     = sprintf('coord_%1$deV_%2$s_%3$s_95_newN2.mat',eVrange,datatype,uncertainty);
savename     = sprintf('coord_%1$deV_%2$s_%3$s_95_thierry.mat',eVrange,datatype,uncertainty);

% Scan settings
start_decade = -1;
stop_decade  = 4;
min_sin2T4   = 0.00001;              % Lower bound for sin2(th4) for the plot

p   = 0.00001;                       % Newton gradient step
err = 0.01;                          % Newton convergence tolerance


% Variables
m4_Y    = [];
sith4_X = [];                       % Plot variables
chi_Z   = [];

% Scanning range
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

% FSD
switch datatype
    case 'Real'
        %sibille = 'SibilleFull';
        sibille = 'Sibille0p5eV';
        chilvl  = chilvl;
    case 'Twin'
        sibille = 'Sibille0p5eV';
end

SysEffects = struct(...
    'RF_EL','OFF',...       % Response Function(RF) EnergyLoss
    'RF_BF','OFF',...       % RF B-Fields
    'RF_RX','OFF',...       % Column Density, inel cross ection
    'FSD','OFF',...         % final states
    'TASR','OFF',...        % tritium activity fluctuations
    'TCoff_RAD','OFF',...   % radiative thoretical corrections (this has to be OFF for KNM1, because included in model)
    'TCoff_OTHER','OFF',... % other theo corr.
    'DOPoff','OFF',...      % doppler effects --> OFF, because in model
    'Stack','OFF',...       % stacking / HV fluctuations
    'FPDeff','OFF');        % detector efficiency

%% Constructor Initialisation
R = MultiRunAnalysis('RunList','KNM1',...
            'chi2',chi2_type,...
            'DataType',datatype,...
            'fixPar','E0 Norm Bkg',...
            'RadiativeFlag','ON',...
            'NonPoissonScaleFactor',NPfactor,...
            'minuitOpt','min ; migrad',...
            'FSDFlag',sibille,...
            'ELossFlag','KatrinT2',...
            'SysBudget',22,...
            'SynchrotronFlag','OFF',...
            'AngularTFFlag','OFF');

% R.ComputeCM('SysEffects',SysEffects,...
%     'BkgCM','ON');

% %% Check systematics errorbars
% times = (R.ModelObj.qUfrac*R.ModelObj.TimeSec);
% 
% CM  = R.FitCM_Obj.CovMatFrac;
% 
% IS = R.ModelObj.TBDIS;
% stat  = sqrt(IS);
% stat  = stat./times;
% 
% err = zeros(1,length(CM));
% for k = (1:length(CM))
%     err(k) = sqrt(CM(k,k));
% end
% 
% prlB = [50 148 216]/255;
% 
% plot(R.ModelObj.qU-18574,err,'color',prlB,'LineWidth',3)
% hold on
% plot(R.ModelObj.qU-18574,stat,'--','color',prlB,'LineWidth',3)
% 
% % Plot style
% xlabel('Retarding energy - 18574 (eV)');
% ylabel('Errorbar');
% legend({'Stat','Syst'},'Location','southwest','box','off');
% 
% % xlim([-90 max(qUc+5)]);
% PRLFormat;
% set(gca, 'YScale', 'log');
% axis square

% R.ModelObj.ComputeTBDDS();
% YD=R.ModelObj.TBDDS;
% R.ModelObj.ComputeTBDIS();
% YI = R.ModelObj.TBDIS;
% R.InitModelObj_Norm_BKG
% R.RunData.TBDIS
% sum(YI)

%% Delta-X2
X0     = fit_chi(0,0,R,eVrange);
chilvl = chilvl + X0;

%% Loop
% Initialisation
%progressbar()
%progressbar(0)
c=0;

diary 'progress.txt'
fprintf('\n==== STARTING SCAN ====\n');
diary off

% Scanning
for m = m4
    % Initializing Newton
    s=1;
    X2_a = fit_chi(m,s,R,eVrange);
    
    diary 'progress.txt'
    fprintf('\nSterile mass : %.3f eV2',m)
    diary off
    
    if (X2_a>chilvl+err)
        diary 'progress.txt'
        fprintf('\n\tReducing \\chi^2 ...\n')
        diary off
        
        % Getting close enough to the value, otherwise minuit is chaotic
        X2_a2=X2_a
        while X2_a2>200
            s=s/1.3;
            X2_a2=fit_chi(m,s,R,eVrange);
        end
        
        diary 'progress.txt'
        fprintf('\tStarting Newton convergence\n')
        diary off
        
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

            % Actualising plot variables
            if ((abs(X2_a-chilvl)<err|abs(X2_b-chilvl)<err) & s>min_sin2T4 & s<=1)
                sith4_X = [sith4_X,s];
                m4_Y    = [m4_Y,m];
                chi_Z   = [chi_Z,X2_a];
            end
        end
    else
        diary 'progress.txt'
        fprintf('\n\tOut of bonds\n')
        diary off
    end
    
    % Progress bar
    c=c+1;
    %progressbar(c/length(m4))
    progress=100*c/length(m4);
    %waitbar(progress,probar,string(progress*100))
    
    diary 'progress.txt'
    Progress = 100*c/length(m4)
    diary off
end

t=toc;
diary 'progress.txt'
fprintf('SCANNING OVER\n\n')
fprintf('Time spent : %.1f hours\n',t/3600)
diary off

%% Datasave
diary 'progress.txt'
fprintf('Saving file ...\n')
diary off

filepath   = [getenv('SamakPath'),'ksn1ana/contour/'];
file       = [filepath,savename];
MakeDir(filepath)
%file_sith4 = '/home/iwsatlas1/guennic/Desktop/Samak2.0/ksn1ana/contour/coord_0to4_dV_90eV_Final';
save(file,'sith4_X','m4_Y','chi_Z');


diary 'progress.txt'
fprintf('==== FINISHED ====\n')
diary off



%% Other
% %% sin(2th)
% 
% sith42_X=[];
% 
% for elt=sith4_X
%     sith42_X = [sith42_X,1-(1-2*elt)^2];
% end
% 
% %% Plots
% 
% %% KATRIN
% plot(sith4_X,m4_Y)
% 
% xlabel('sin^2(\theta_4)');
% ylabel('m_4^2');
% legend({'X^2 = 4.6'},'Location','southwest')
% 
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% 
% %% GIUNTI
% figure;
% [ms,cut]=max(sith42_X); N=length(sith42_X);
% %scatter(sith42_X(cut:N),m4_Y(cut:N))      % Remove the bottom tail
% plot(sith42_X,m4_Y)
% 
% xlabel('sin^2(2\theta_4)');
% ylabel('m_4^2');
% legend({'X^2 = 4.6'},'Location','southwest')
% 
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
