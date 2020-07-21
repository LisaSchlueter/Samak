%
% Th. Lasserre, 2014
% Class To Compute Relic Neutrino Capture on Tritium
% Rate and Spectrum in KATRIN
%
% CEA - Saclay
%

classdef TritiumRelicNu

    properties (Constant)
        sin2thetaW = 0.23160;     % Weinberg Angle PDG2010
        alphaem = 1/137.035999074;     % Fine Structure Cte
        me = 510.998918;         % Electron Mass keV
        c = 299792458;           % Speed of Light m.s^-1
        hbarc = 197.3269631e-15; % Planck Reduced Cte MeV.m
        Gf = 1.166637e-11*...
            (197.3269631e-15)^3; %Fermi Cte  Mev.m^3 PDG2010
        Mass_n = 939565.360;     % Neutron mass keV
        Mass_p = 938272.029;     % Proton mass keV
        MH3   = 3.0160492 * 931494.061; % Tritium mass keV 3.0160492787
        MHe3  = 3.0160293 * 931494.061; % Helium-3 mass keV 3.0160293217
        Gfbis  = 1.1663787e-11;       % Fermi constant Mev-2, PDG2012
        coscab = 0.9743;              % cossqrcabibbo angle, PDG 2012
        gA     = 1.247; gV     = 1.000;
        Avogadro = 6.02214129e23 % mol?1
    end
   
    properties (Access=public)
        % Detector Properties
        KatrinDetEff=0.9;          % Detector Efficiency
        ThetaMax=50.77;            % Maximum Acceptance Angle (deg)
        CosThetaMax=0.6324;        % Cos(ThetaMax)
        Tpurity=0.95;              % Tritium Purity in Katrin
        Rs=4.1;                    % Radius Flux Tube in cm
        Rhod=5e17;                 % Column Density in molecules/cm^2
        SourceType;                % 1 for Tritium
        WGTSbetaDecay=9.42e10;     % Number of Beta Decay per sec in WGTS
        
        % Tritium Decay
        Z=2;                       % 3He
        Q=18.574;                  % Tritium Q-value, keV, Regular
        tau=12.32*86400*365.25;    % Tritium Tau in Sec
        TritiumDecayConst=1.78283e-09; % s^1
        TritiumActivity_perMol  = 1.0736e+15; % Decay per Sec per Mol
        TritiumActivity_perGram = 3.5598e+14; % Decay Per Sec Per Gram
    end
    
    properties (Dependent=true,Access=public)
    end
    
    properties (Access=public)
        Eta;     % relic neutrino overdensity
        tex;     % Experiment time in s
        Ee0;     % End point expressed in total neutrino energy
        nE;      % number of energy bins
        nF;      % number of fine energy bins
        mnu=0;   % neutrino 1, 2, 3 masses
        sigma_E; % Energy resolution
        Temin;   % Minimum e kinetic energy
        Temax;   % Maximum e kinetic energy
        nTe;     % Kinetic energy Bins
        TeStep;  % energy bin steps width
        A;       % Number of Nucleons
        F;       % Katrin Tritium Beta Decay Parameter, Fermi Function
        R;       % Katrin Tritium Beta Decay Parameter
        mate;    % energy resolution
        
        % Relic Neutrinos Spectrum
        RNS_nTnu;          % Relic Neutrino Spectrum bin width (keV)
        TnuStep;           % energy bin steps width
        RNS_Pnu ;          % Cosmic neutrino momentum
        RNS_Enu; RNS_EnuMin ; RNS_EnuMax ; % Cosmic neutrino energy
        RNS_EnuError;      % neutrino capture : neutrino energy Error
        RNS_EnuFluct;      % neutrino capture : neutrino energy fluct
        RNS_EnuFluctError; % neutrino capture : neutrino energy fluct error
        NormRNS = 56;      % Nbx neutrino per specie per cm3, 112 including antiparticle
        RNS_Snu;           % Cosmic neutrino spectrum
        RNS_Tnu = 1.67e-7; % RNS Temperature (1.9 K - 1.67e-7 KeV)
        
        % Relic Neutrinos: Electron Capture
        RNS_Pe ;     % neutrino capture : electron momentum
        RNS_Ee;      % neutrino capture : electron energy
        RNS_EeError; % neutrino capture : electron energy Error
        RNS_beta ;   % neutrino capture : beta factor
        RNS_F;       % neutrino capture : Fermi Function
        RNS_EeFluct; % neutrino capture : electron energy fluct
        RNS_EeFluctError; % neutrino capture : electron energy fluct error
        SigmaV;      % cross Section times Velocity
        Tmass;       % Tritium Mass in Gram
        Tmol;        % Tritium Number of Mole
        RateCaptureT;% capture rate on tritiulm per mol per year
        RateCaptureT_KATRIN; % Capture rate in Katrin
        RateCaptureT_KATRINerror; % its uncertainty, poissonian
        T_eCapture_S;% Electron Capture Spectrum on Tritium
        Te;          % Electron kinetic energy
        
        % Fit
        TeFirstBin; % First Bin for Fit purpose
        TeLastBin;  % Last Bin for Fit purpose
                       
        % Fermi Function
        FermiFuncType; % 0 Non Relativistic , 1 Relativistic
        FermiF;
        
        % Normalization Type
        NormType;     % Type of Normalization
        TotalCounts;  % Given Total Number of Counts   

    end
     
     methods
         function obj = TritiumRelicNu(varargin)
             obj.nF = 100; % 100 eV bins
             p = inputParser;
             p.addParameter('NormType','TMASS',@(x)ismember(x,{'TMASS','COUNTS'}));
             p.addParameter('Tmass',40,@(x)isfloat(x) && x>0);
             p.addParameter('TotalCounts',1e9,@(x)isfloat(x) && x>-1);
             p.addParameter('Eta',1,@(x)isfloat(x));
             p.addParameter('activity',1,@(x)isfloat(x) && x>0);
             p.addParameter('tex',1,@(x)isfloat(x) && x>0);
             p.addParameter('RNS_EnuMin',-9,@(x)isfloat(x) && x>0);
             p.addParameter('RNS_EnuMax',-5,@(x)isfloat(x) && x>0);
             p.addParameter('RNS_nTnu',100,@(x)isfloat(x) && x>0);
             p.addParameter('Temin',18.5,@(x)isfloat(x) && x>0);
             p.addParameter('Temax',18.6,@(x)isfloat(x) && x>0);
             p.addParameter('nTe',100,@(x)isfloat(x) && x>0);
             p.addParameter('mnu',0,@(x)isfloat(x));
             p.addParameter('TeFirstBin',1,@(x)isfloat(x));
             p.addParameter('TeLastBin',100,@(x)isfloat(x));
             p.addParameter('energy_resol',1e-99,@(x)isfloat(x) && x>=0);
             
             p.parse(varargin{:});
             obj.NormType=p.Results.NormType;
             obj.Tmass=p.Results.Tmass;
             obj.TotalCounts=p.Results.TotalCounts;
             obj.Eta = p.Results.Eta;
             obj.tex=p.Results.tex;% *365.25*86400;
             obj.mnu=p.Results.mnu;
             
             obj.RNS_EnuMin = p.Results.RNS_EnuMin;
             obj.RNS_EnuMax = p.Results.RNS_EnuMax;
             obj.RNS_nTnu=p.Results.RNS_nTnu;
             obj.RNS_Enu=logspace(obj.RNS_EnuMin,obj.RNS_EnuMax,obj.RNS_nTnu)';
             obj.TnuStep = obj.RNS_Enu(2)-obj.RNS_Enu(1);
             
             obj.Temin = p.Results.Temin;
             obj.Temax = p.Results.Temax;
             obj.nTe=p.Results.nTe;
             obj.Te=linspace(obj.Temin,obj.Temax,obj.nTe)';
             obj.TeStep = obj.Te(2)-obj.Te(1);
             obj.sigma_E = p.Results.energy_resol;

             obj.TeFirstBin=p.Results.TeFirstBin;
             obj.TeLastBin=p.Results.TeLastBin;
            
             % Energy Bin Step
             obj.TnuStep = obj.RNS_Enu(2)-obj.RNS_Enu(1);
             
             % Define Energy Resolution
             obj.mate=eres(obj.RNS_Enu,obj.RNS_EnuMin,obj.RNS_EnuMax,obj.RNS_nTnu,obj.sigma_E);
                                                                     
             % ----------------------------------------------------------
             %  Relic Neutrino Spectrum
             % ----------------------------------------------------------
             obj.RNS_Pnu = obj.RNS_Enu;
             obj.RNS_Snu = obj.RNS_Pnu.^2.*(exp(obj.RNS_Pnu./(obj.RNS_Tnu))+1).^-1;
             obj.RNS_Snu = obj.RNS_Snu./sum(obj.RNS_Snu)*obj.NormRNS; % - cm^-3 keV^-1
             
             % ----------------------------------------------------------
             % Kinematic Variable Neutrino Capture
             % ----------------------------------------------------------
             obj.RNS_Ee  = obj.Q + obj.me + obj.RNS_Enu.*obj.RNS_Snu;
             obj.RNS_Pe = real(obj.RNS_Ee.^2-obj.me.^2).^.5;  
             obj.RNS_beta = obj.RNS_Pe./obj.RNS_Ee;  
             
             % ----------------------------------------------------------
             %  Fermi Function
             %  Case 0 : Non Relativistic J.J. Simpson, Phys. Rev. D23 (1981) 649
             % ----------------------------------------------------------
             obj.A=2*pi.*obj.alphaem.*obj.Z./obj.RNS_beta;
             obj.RNS_F=obj.A./(1-exp(-obj.A));
             
             % Cross Section
             C = 0.9987^2 + 1.2695 * sqrt(3) * 0.964 ;
             obj.SigmaV = obj.Gf ./ pi .* obj.RNS_Pe .* obj.RNS_Ee .* obj.RNS_F * C;

             % Nbx of Interaction in Tritium
             switch obj.NormType
                 case 'COUNTS'
                     % ----------------------------------------------------------
                     % Normalisation based on a given total number of counts
                     % ----------------------------------------------------------
                     obj.NormRNS = obj.TotalCounts;
                     obj.RateCaptureT_KATRIN = obj.TotalCounts;
                 case 'TMASS'
                     obj.NormRNS = 1;
                     obj.RateCaptureT = obj.Eta .* 2.85e-2 .* obj.SigmaV ./ 1e-45; % per year per mol
                     obj.Tmol = obj.Tmass / 3.0160492;
                     obj.RateCaptureT_KATRIN = obj.RateCaptureT .* obj.Tmol .* obj.tex;
             end
             
             % Neutrino Capture: Electron Spectrum
             % ----------------------------------------------------------
             % Apply Energy Resolution
             % Integrate over the bin step size (approximation)
             % ----------------------------------------------------------
             obj.T_eCapture_S = zeros(obj.nTe,1);
             if obj.sigma_E < 1e-9
                 clear k; k=find(obj.Te>(obj.Q+obj.mnu-1e-5) & obj.Te<(obj.Q+obj.mnu+1e-5));
                 obj.T_eCapture_S(k) = obj.RateCaptureT_KATRIN(1);
             end
             if obj.sigma_E > 1e-9
                 clear tmp ; tmp = normpdf(obj.Te,obj.Q+obj.mnu,obj.sigma_E);
                 obj.T_eCapture_S = tmp/sum(tmp)*obj.RateCaptureT_KATRIN(1);
                 %obj.T_eCapture_S = obj.mate * obj.T_eCapture_S;
             end
             
             % Spectrum Errors
             obj.RateCaptureT_KATRINerror = real(obj.RateCaptureT_KATRINerror).^.5;
       
             % Statistical Fluctuations & Errors
%              obj.RNS_EeFluct = ...
%                   obj.RNS_Ee + ...
%                   obj.RNS_EeError.*randn(obj.RNS_nTnu,1);
%              obj.RNS_EeFluctError  = obj.RNS_EeFluct.^0.5;
                                       
         end % constructor
         
     end % methods
     
     methods
             % methods here
     end 

end % class

function r = eres(E,Emin,Emax,Nbins,sigma)
% Produces Energy resolution matrix from Ebins to Nbins:
% r is [nB nE].
Estep = (Emax-Emin)/Nbins;
dn = repmat(reshape(linspace(Emin,Emax-Estep,Nbins),[Nbins 1]),[1 numel(E)]);
up = repmat(reshape(linspace(Emin+Estep,Emax,Nbins),[Nbins 1]),[1 numel(E)]);
f=@(x,a,b,sigma).5*double((erf((b-x)./(sqrt(2)*sigma))-erf((a-x)./(sqrt(2)*sigma))));
r = f(repmat(reshape(E,[1 numel(E)]),[Nbins 1]),dn,up,sigma);
end