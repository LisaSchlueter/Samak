addpath(genpath('../../../Samak2.0'));

MultiKrObj = InitKrKATRIN_Lisa();
MultiKrObj.ComputeKrDS;
figure(1);
plot(MultiKrObj.Te*1e-3, MultiKrObj.KrDS);
ylabel('$$\frac{\textrm{dN}}{\textrm{dE}}$$','Interpreter','latex');
xlabel('$$\textrm{E}_{\textrm{kin}} \textrm{[keV]}$$', 'Interpreter', 'latex');
xlim([30.45 30.480]);
grid on;
set(gcf,'Color','w') %background white
export_fig ../krypton_L3Satellites_Lisa/plots/SatellitePeaksDS_L3_32_new.pdf
%figure(3); %time distribution (qU and qUfrac from KATRIN_Lisa.m)  
%bar(MultiKrObj.qU*1e-3, MultiKrObj.qUfrac);

figure(2);
MultiKrObj.ComputeKrIS;
plot(MultiKrObj.qU*1e-3 , MultiKrObj.KrIS);
grid on;
xlabel('retarded potential qU [keV]');
ylabel('N');
xlim([30.45 30.480]);
export_fig ../krypton_L3Satellites_Lisa/plots/SatellitePeaksIS_L3_32_new.pdf

