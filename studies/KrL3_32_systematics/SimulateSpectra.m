function SamakSpec = SimulateSpectra(varargin)
 p = inputParser;
 p.addParameter('KrObject','',@(x) isa(x,'Kr'));
 p.addParameter('StatFluct', 'ON', @(x)ismember(x,{'ON', 'OFF'}));
 p.parse(varargin{:});
 obj = p.Results.KrObject;
StatFluct = p.Results.StatFluct;

addpath(genpath('../../../Samak2.0'));

%obj= InitKrKATRIN_krl332();
obj.L3_32_W_i =1.557; %approx value of L3-32 Data
obj.L3_32_E_i = 30475.0;
obj.TimeSec =10*60*60; %

obj.ComputeKrDS;
obj.ComputeKrIS;
%if strcmp ('StatFluct', 'ON') 
%obj.AddStatFluctKrIS(); % end %my true parameters + stat fluctuations

SamakSpec = [flip(obj.qU) , flip(obj.KrIS) , flip(obj.KrISE)];
%save

%errorbar(obj.qU, obj.KrIS, obj.KrISE,'x')
filepath = sprintf('../../krypton-data/33051_L3-32_Simulation_nofluct/33051_channel%u.txt',(obj.FPD_Pixel-1));
save(filepath, 'SamakSpec', '-ascii');
end
