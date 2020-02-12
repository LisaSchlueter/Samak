%
% Read / Load Magnetic Field as a function of z
%

Bz=load('magnetic_field_on_particle_axis_v2.mat','z','B');
z=Bz.z;
B=Bz.B;

figure(1)
plt1= Plot(z,B);
plt1.LineWidth = 2;
plt1.LineStyle = {'-'};
plt1.Markers   = {'.'};
plt1.MarkerSpacing = 1;
pltstr         = sprintf('');
plt1.Title     = pltstr; 
plt1.XLabel    = 'z (m)';  
plt1.YLabel    = 'B (T)';  
plt1.YScale    = 'lin';  
plt1.XScale    = 'lin';  
plt1.FontSize  = 16;
plt1.Legend    = {'B field along z axis'};
plt1.LegendLoc = 'NorthWest';
plttitle       = sprintf('KATRIN_BversuZ.png');
plt1.export(plttitle);
