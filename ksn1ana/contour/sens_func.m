function sens_func(contour_settings)
    %% KATRIN Sterile Neutrino sensitivity - KSN1
    %  Nathan Le Guennic - 2020

    tic;
    %% Settings
    % Parameters

    CL           = contour_settings.CL;                  % Confidence level                 % 90% - 95% - 99%
    datatype     = contour_settings.datatype;                                               % Real - Twin
    uncertainty  = contour_settings.uncertainty;                                            % syst - stat
    NPfactor     = contour_settings.NPfactor;           % Non poisson factor                % 1 - 1.064
    d            = contour_settings.scan_step;          % Number of dots per decade
    eVrange      = contour_settings.eVrange;            % eV below the endpoint
    ActiveNeut   = contour_settings.activeFlag;         % Activate the active neutrino fit  % FREE OFF FIX

    % Name for the datafile
    savename     = sprintf('coord_%1$deV_%2$s_%3$s_%4$dCL',eVrange,datatype,uncertainty,CL);

    
    %% Scan settings
    start_decade = -1;
    stop_decade  = 4;
    min_sin2T4   = 0.00001;                 % Lower bound for sin2(th4) for the plot

    p   = 0.00001;                          % Newton gradient step
    err = 0.01;                             % Newton convergence tolerance


    % Variables
    m4_Y    = [];
    sith4_X = [];                           % Plot variables
    chi_Z   = [];
    m_beta  = [];

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
            sibille = 'SibilleFull';
            chilvl  = chilvl;
        case 'Twin'
            sibille = 'Sibille0p5eV';
    end

    % Active neutrino
    switch ActiveNeut
        case 'FREE'
            free_para = 'mNu E0 Norm Bkg';
            savename  = [savename,'_freeM.mat'];
        case 'OFF'
            free_para = 'E0 Norm Bkg';
            savename  = [savename,'.mat'];
        case 'FIX'
            free_para = 'E0 Norm Bkg';
            savename  = [savename,'_fixM.mat'];
    end


    %% Constructor Initialisation
    R = MultiRunAnalysis('RunList','KNM1',...
                'chi2',chi2_type,...
                'DataType',datatype,...
                'fixPar',free_para,...
                'RadiativeFlag','ON',...
                'NonPoissonScaleFactor',NPfactor,...
                'minuitOpt','min ; migrad',...
                'FSDFlag',sibille,...
                'ELossFlag','KatrinT2',...
                'SysBudget',22);
            
    %% Delta-X2
    X0     = fit_chi(0,0,R,eVrange);
    chilvl = chilvl + X0;
%     X0   = 0;

    %% Loop
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

                    if strcmp(ActiveNeut,'FREE')
                        m_beta  = [m_beta,R.FitResult.par(1)]; % When m_beta is activated
                    end

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
    
    %% Convert data
    % Delta M
    switch ActiveNeut
        case 'FREE'
            DM2 = m4_Y-m_beta;
        case 'OFF'
            DM2 = m4_Y;
        case 'FIX'
            DM2 = m4_Y+0.956847;    % Best KN1 fit
    end
    
    % sin(theta) - sin(2 theta)
    si2th4_X = 1-(1-2*sith4_X).^2;
    
    %% Datasave
    diary 'progress.txt'
    fprintf('Saving file ...\n')
    diary off

    filepath   = [getenv('SamakPath'),'ksn1ana/contour/contour_files/V2/'];
    file       = [filepath,savename];
    MakeDir(filepath)

    save(file,'contour_settings','R','sith4_X','si2th4_X','m4_Y','m_beta','DM2','chi_Z','X0');

    diary 'progress.txt'
    fprintf('==== FINISHED ====\n')
    diary off
end