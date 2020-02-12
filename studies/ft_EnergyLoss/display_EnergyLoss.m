function display_EnergyLoss(varargin)% Display E-Loss Functions
% T. Lasserre, 14/12/2018
Nscatt = 8;

% Parser
p = inputParser;
p.addParameter('Nscatt',Nscatt,@(x)isfloat(x) && x>0);
p.parse(varargin{:});
Nscatt                 = p.Results.Nscatt;

% Create WGTS objects
Aseev         = WGTSMACE('ELossFlag','Aseev');
Abdurashitov  = WGTSMACE('ELossFlag','Abdurashitov');
CW_GLT        = WGTSMACE('ELossFlag','CW_GLT');
CW_G2LT       = WGTSMACE('ELossFlag','CW_G2LT');
% Init
Aseev.InitializeELossFunction;
Abdurashitov.InitializeELossFunction
CW_GLT.InitializeELossFunction
CW_G2LT.InitializeELossFunction
% Get cells of eloss function handles 
cells_Aseev        = Aseev.ComputeELossFunction();
cells_Abdurashitov = Abdurashitov.ComputeELossFunction();
cells_CW_GLT       = CW_GLT.ComputeELossFunction();
cells_CW_G2LT      = CW_G2LT.ComputeELossFunction();
% Get 1st inelastic scattering eloss functions
Aseev1        = cells_Aseev{Nscatt};
Abdurashitov1 = cells_Abdurashitov{Nscatt};
CW_GLT1       = cells_CW_GLT{Nscatt};
CW_G2LT1      = cells_CW_G2LT{Nscatt};
%% Plot 1st inelastic scattering eloss functions 
e=-1000:.1:1000;
pel = Plot(e,Aseev1(e),e,Abdurashitov1(e),e,CW_GLT1(e),e,CW_G2LT1(e));
pel.YGrid = 'on'; % 'on' or 'off'
pel.BoxDim = [14, 6]; %[width, height] in inches
pel.Legend = {'Aseev','Abdurashitov','CW_GLT','CW_G2LT'}; % legends,
pel.XLabel = 'Energy Loss(eV)'; % xlabel
pel.YLabel = 'pdf'; %ylabel
pel.Title  = sprintf('Energy Loss Candidate Functions - Inelastic Scattering N° %0.f',Nscatt); % plot title
pel.LineWidth = [2, 2, 2, 2]; % three line widths
pel.LineStyle = {'-', '-', '-', '-'}; % three line styles
pel.MarkerSpacing = [10, 10 , 10 , 10];
pel.YScale = 'lin'; % 'linear' or 'log'
pel.XLim   = [0 50*Nscatt];