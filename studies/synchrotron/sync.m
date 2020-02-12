%
% Compute Synchrotron Energy Loss
%
close all;

ME      = 510998.0;           % eV
CLIGHT  = 299792458;                % m per second
BSOURCE = (2.52);             % Bfield at source in Tesla
ZREAR   = -44.0;              % Rear of WGTS
ZSTART  = -39.0;              % Middle of WGTS
ZFRONT  = -34.0;              % Front of WGTS
ZEND    = -15.0;              % Spectrometer AP
E_source = 18575.0;           % Electron energy keV

% Read B field
Bz=load('magnetic_field_on_particle_axis_v2.mat','z','B');
global LookUpz; LookUpz = Bz.z';
global LookUpB; LookUpB = Bz.B';

%for i=0:1:204
for i=0:1:210
         theta_source(i+1) = i*0.25/180.0*pi;   
         g = @(x)tof_integrand(x,E_source,BSOURCE,theta_source(i+1));
         DeltaERear(i+1)   = integral(g,ZREAR,ZEND);
         DeltaEMiddle(i+1) = integral(g,ZSTART,ZEND);
         DeltaEFront(i+1)  = integral(g,ZFRONT,ZEND);
fprintf('angle\t%.3f degree \t synchrotron loss %g \t meV\n',(theta_source(i+1)/pi*180.0), DeltaEMiddle(i+1)*1e3);
end

%% Averaged Energy Loss
MeanDeltaERear     = mean(DeltaERear)*1e3;
MeanDeltaEMiddle   = mean(DeltaEMiddle)*1e3;
MeanDeltaEFront    = mean(DeltaEFront)*1e3;   
fprintf('Averaged Rear: %g meV\n',MeanDeltaERear);
fprintf('Averaged Middle: %g meV\n',MeanDeltaEMiddle);
fprintf('Averaged Front: %g meV\n',MeanDeltaEFront);

%%
figure(1)
plt1= Plot((theta_source/pi*180.0),DeltaERear*1e3,...
(theta_source/pi*180.0),DeltaEMiddle*1e3,...
(theta_source/pi*180.0),DeltaEFront*1e3);...
plt1.LineWidth = 2;
plt1.LineStyle = {'-','-','-'};
plt1.Markers   = {'.','.','.'};
plt1.MarkerSpacing = 1;
pltstr         = sprintf('');
plt1.Title     = pltstr; 
plt1.XLabel    = 'pitch angle (degree)';  
plt1.YLabel    = 'synchrotron loss (meV)';  
plt1.YScale    = 'lin';  
plt1.XScale    = 'lin';  
plt1.FontSize  = 16;
plt1.Legend    = {'Start in Rear of WGTS','Start in Middle of WGTS','Start in Front of WGTS'};
plt1.LegendLoc = 'NorthWest';
plttitle       = sprintf('Synchrotron-WGTS.png');
plt1.export(plttitle);

figure(2)
plt1= Plot((theta_source/pi*180.0),DeltaERear./DeltaEMiddle,...
(theta_source/pi*180.0),DeltaEMiddle./DeltaEMiddle,...
(theta_source/pi*180.0),DeltaEFront./DeltaEMiddle);...
plt1.LineWidth = 2;
plt1.LineStyle = {'-','-','-'};
plt1.Markers   = {'.','.','.'};
plt1.MarkerSpacing = 1;
pltstr         = sprintf('');
plt1.Title     = pltstr; 
plt1.XLabel    = 'pitch angle (degree)';  
plt1.YLabel    = 'synchrotron loss (meV)';  
plt1.YScale    = 'lin';  
plt1.XScale    = 'lin';  
plt1.FontSize  = 16;
plt1.Legend    = {'Start in Rear of WGTS','Start in Middle of WGTS','Start in Front of WGTS'};
plt1.LegendLoc = 'NorthWest';
plttitle       = sprintf('Synchrotron-WGTS.png');
plt1.export(plttitle);

%
%%
figure(3)
plt1= Plot(sin(theta_source),DeltaERear*1e3,...
(sin(theta_source)),DeltaEMiddle*1e3,...
(sin(theta_source)),DeltaEFront*1e3);...
plt1.LineWidth = 2;
plt1.LineStyle = {'-','-','-'};
plt1.Markers   = {'.','.','.'};
plt1.MarkerSpacing = 1;
pltstr         = sprintf('');
plt1.Title     = pltstr; 
plt1.XLabel    = 'pitch angle (degree)';  
plt1.YLabel    = 'synchrotron loss (meV)';  
plt1.YScale    = 'lin';  
plt1.XScale    = 'lin';  
plt1.FontSize  = 16;
plt1.Legend    = {'Start in Rear of WGTS','Start in Middle of WGTS','Start in Front of WGTS'};
plt1.LegendLoc = 'NorthWest';
plttitle       = sprintf('Synchrotron-WGTS.png');
plt1.export(plttitle);