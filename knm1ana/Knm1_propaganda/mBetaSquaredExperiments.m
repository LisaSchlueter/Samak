function [mBetaSquared , mBetaSquaredPDG] = mBetaSquaredExperiments()

%% PDG Average
mBetaSquaredPDG = struct('Label',{},'Year',{},'Reference',{},'mBetaSquared',{},'Tot',{},'Marker',{},'Color',{});
mBetaSquaredPDG(1).Label='PDG 2019 Average';
mBetaSquaredPDG(1).Year=2019;
mBetaSquaredPDG(1).Reference='Phys. Rev. D 98, 030001 (2018)';
mBetaSquaredPDG(1).mBetaSquared=-0.6;
mBetaSquaredPDG(1).Tot=1.9;
mBetaSquaredPDG(1).Marker='.';
mBetaSquaredPDG(1).Color=rgb('CadetBlue');

%% From PDG2019 - DATA m2
mBetaSquared = struct('Experiment',{},'Label',{},'Year',{},'Reference',{},'mBetaSquared',{},'Stat',{},'Sys',{},'Tot',{},'Marker',{},'Color',{});

mBetaSquared(1).Experiment='Los Alamos';
mBetaSquared(1).Label='Robertson 1991';
mBetaSquared(1).Year=1990.5;
mBetaSquared(1).Reference='Phys.Rev.Lett. 67 (1991) 957-960';
mBetaSquared(1).mBetaSquared=-147;
mBetaSquared(1).Stat=68;
mBetaSquared(1).Sys=41;
mBetaSquared(1).Tot=sqrt(mBetaSquared(1).Stat^2+mBetaSquared(1).Sys^2);
mBetaSquared(1).Marker='+';
mBetaSquared(1).Color=rgb('DarkSlateGray');

mBetaSquared(2).Experiment='Tokyo';
mBetaSquared(2).Label='KAWAKAMI 1991';
mBetaSquared(2).Year=1991.5;
mBetaSquared(2).Reference='Phys.Lett. B256 (1991) 105-111';
mBetaSquared(2).mBetaSquared=-65;
mBetaSquared(2).Stat=85;
mBetaSquared(2).Sys=65;
mBetaSquared(2).Tot=sqrt(mBetaSquared(2).Stat^2+mBetaSquared(2).Sys^2);
mBetaSquared(2).Marker='<';
mBetaSquared(2).Color=rgb('Indigo');

mBetaSquared(3).Experiment='Zurich';
mBetaSquared(3).Label='Holzschuh 1992';
mBetaSquared(3).Year=1992;
mBetaSquared(3).Reference='Phys.Lett. B287 (1992) 381-388';
mBetaSquared(3).mBetaSquared=-24;
mBetaSquared(3).Stat=48;
mBetaSquared(3).Sys=61;
mBetaSquared(3).Tot=sqrt(mBetaSquared(3).Stat^2+mBetaSquared(3).Sys^2);
mBetaSquared(3).Marker='*';
mBetaSquared(3).Color=rgb('Brown');

mBetaSquared(4).Experiment='Mainz';
mBetaSquared(4).Label='Weinheimer 1993';
mBetaSquared(4).Year=1992.5;
mBetaSquared(4).Reference='Phys.Lett. B300 (1993) 210-216';
mBetaSquared(4).mBetaSquared=-39;
mBetaSquared(4).Stat=34;
mBetaSquared(4).Sys=15;
mBetaSquared(4).Tot=sqrt(mBetaSquared(4).Stat^2+mBetaSquared(4).Sys^2);
mBetaSquared(4).Marker='d';
mBetaSquared(4).Color=rgb('CadetBlue');

mBetaSquared(5).Experiment='China';
mBetaSquared(5).Label='Sun 1993';
mBetaSquared(5).Year=1993.5;
mBetaSquared(5).Reference='Chin.J.Nucl.Phys. 15 (1993) 261';
mBetaSquared(5).mBetaSquared=-31;
mBetaSquared(5).Stat=75;
mBetaSquared(5).Sys=48;
mBetaSquared(5).Tot=sqrt(mBetaSquared(5).Stat^2+mBetaSquared(5).Sys^2);
mBetaSquared(5).Marker='v';
mBetaSquared(5).Color=rgb('Gold');

mBetaSquared(6).Experiment='LLNL';
mBetaSquared(6).Label='Stoeffl 1995';
mBetaSquared(6).Year=1994.5;
mBetaSquared(6).Reference='Phys.Rev.Lett. 75 (1995) 3237-3240';
mBetaSquared(6).mBetaSquared=-130;
mBetaSquared(6).Stat=20;
mBetaSquared(6).Sys=15;
mBetaSquared(6).Tot=sqrt(mBetaSquared(6).Stat^2+mBetaSquared(6).Sys^2);
mBetaSquared(6).Marker='>';
mBetaSquared(6).Color=rgb('Brown');

mBetaSquared(7).Experiment='Troitsk';
mBetaSquared(7).Label='Belesev 1995';
mBetaSquared(7).Year=1995.5;
mBetaSquared(7).Reference='Phys.Lett. B350 (1995) 263-272';
mBetaSquared(7).mBetaSquared=-22;
mBetaSquared(7).Stat=sqrt(4.8^2/2);
mBetaSquared(7).Sys=sqrt(4.8^2/2);
mBetaSquared(7).Tot=4.8;
mBetaSquared(7).Marker='s';
mBetaSquared(7).Color=rgb('IndianRed');

mBetaSquared(8).Experiment='Mainz';
mBetaSquared(8).Label='Weinheimer 1999';
mBetaSquared(8).Year=1998.5;
mBetaSquared(8).Reference='Phys.Lett. B460 (1999) 219-226';
mBetaSquared(8).mBetaSquared=-3.7;
mBetaSquared(8).Stat=5.3;
mBetaSquared(8).Sys=2.1;
mBetaSquared(8).Tot=sqrt(mBetaSquared(8).Stat^2+mBetaSquared(8).Sys^2);
mBetaSquared(8).Marker='d';
mBetaSquared(8).Color=rgb('CadetBlue');

mBetaSquared(9).Experiment='Troitsk';
mBetaSquared(9).Label='Lobashev 1999';
mBetaSquared(9).Year=1999.5;
mBetaSquared(9).Reference='Phys.Lett. B460 (1999) 227-235';
mBetaSquared(9).mBetaSquared=-1.9;
mBetaSquared(9).Stat=3.4;
mBetaSquared(9).Sys=2.2;
mBetaSquared(9).Tot=sqrt(mBetaSquared(9).Stat^2+mBetaSquared(9).Sys^2);
mBetaSquared(9).Marker='s';
mBetaSquared(9).Color=rgb('IndianRed');

mBetaSquared(10).Experiment='Mainz';
mBetaSquared(10).Label='Kraus 2005';
%mBetaSquared(10).Year=2000;
mBetaSquared(10).Year=2005;
mBetaSquared(10).Reference='Eur.Phys.J. C40 (2005) 447-468';
mBetaSquared(10).mBetaSquared=-0.6;
mBetaSquared(10).Stat=2.2;
mBetaSquared(10).Sys=2.1;
mBetaSquared(10).Tot=sqrt(mBetaSquared(10).Stat^2+mBetaSquared(10).Sys^2);
mBetaSquared(10).Marker='d';
mBetaSquared(10).Color=rgb('CadetBlue');

mBetaSquared(11).Experiment='Troitsk';
mBetaSquared(11).Label='Aseev 2011';
%mBetaSquared(11).Year=2000.5;
mBetaSquared(11).Year=2011;
mBetaSquared(11).Reference='Phys.Rev. D84 (2011) 112003';
mBetaSquared(11).mBetaSquared=-0.67;
mBetaSquared(11).Stat=1.89;
mBetaSquared(11).Sys=1.68;
mBetaSquared(11).Tot=sqrt(mBetaSquared(11).Stat^2+mBetaSquared(11).Sys^2);
mBetaSquared(11).Marker='s';
mBetaSquared(11).Color=rgb('IndianRed');

% Samak DATA
% mBetaSquared(12).Experiment='KATRIN 1^{st} run: 22% Tritium activity - 23 days ';
% mBetaSquared(12).Label='KATRIN KNM1 2019';
% mBetaSquared(12).Year=2019+6/12;
% mBetaSquared(12).Reference='';
% mBetaSquared(12).mBetaSquared=-0.98;
% mBetaSquared(12).Stat=0.935;  % average+/-
% mBetaSquared(12).Sys=0.2764;  % subtracted in quadrature
% mBetaSquared(12).Tot=sqrt(mBetaSquared(12).Stat^2+mBetaSquared(12).Sys^2);
% mBetaSquared(12).Marker='o';
% mBetaSquared(12).Color=rgb('DarkGreen');

%Samak MC
% mBetaSquared(12).Experiment='KATRIN 1^{st} Science Run Sensitivity (MC)';
% mBetaSquared(12).Label='KATRIN KNM1 2019 Sensitivity (MC)';
% mBetaSquared(12).Year=2019+6/12;
% mBetaSquared(12).Reference='';
% mBetaSquared(12).mBetaSquared=0;
% mBetaSquared(12).Stat=0.935;  % average+/-
% mBetaSquared(12).Sys=0.3;  % subtracted in quadrature
% mBetaSquared(12).Tot=0.79;
% mBetaSquared(12).Marker='o';
% mBetaSquared(12).Color=rgb('DarkGreen');%


% KATRIN FINAL
mBetaSquared(12).Experiment='KATRIN 1^{st} Science Run ';
mBetaSquared(12).Label='KATRIN KNM1 2019';
mBetaSquared(12).Year=2019+6/12;
mBetaSquared(12).Reference='';
mBetaSquared(12).mBetaSquared=-1;
mBetaSquared(12).Stat=0.9474;  % average+/-
mBetaSquared(12).Sys=0.32;  % subtracted in quadrature
mBetaSquared(12).Tot=sqrt(mBetaSquared(12).Stat^2+mBetaSquared(12).Sys^2);
mBetaSquared(12).Marker='o';
mBetaSquared(12).Color=rgb('DarkGreen');

mBetaSquared(13).Experiment='KATRIN sensitivity - 100% T-activity - 5 y ';
mBetaSquared(13).Label='KATRIN 2024';
mBetaSquared(13).Year=2024+6/12;
mBetaSquared(13).Reference='';
mBetaSquared(13).mBetaSquared=0;
mBetaSquared(13).Stat=sqrt(0.04^2/2);
mBetaSquared(13).Sys=sqrt(0.04^2/2);
mBetaSquared(13).Tot=0.04;
mBetaSquared(13).Marker='o';
mBetaSquared(13).Color=rgb('DarkGreen');
