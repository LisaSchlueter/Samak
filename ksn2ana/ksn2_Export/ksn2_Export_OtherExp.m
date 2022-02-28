
% import contours and convert to hdf5 file...
fdir = [getenv('SamakPath'),'SterileAnalysis/GridSearchFiles/Knm1/Others/'];
savename =  sprintf('%sContours_OtherExp.h5', fdir);
system(sprintf('rm %s',savename));

%% RAA
f_RAAa = importdata([fdir,'RAA2011_95a.txt']);
f_RAAb = importdata([fdir,'RAA2011_95b.txt']);
f_RAAc = importdata([fdir,'RAA2011_95c.txt']);

ref = "https://arxiv.org/abs/1101.2755";
cl = 95;

h5create(savename,'/RAA/Reference',size(ref),'Datatype','string');
h5write(savename,'/RAA/Reference',ref)
h5create(savename,'/RAA/CL',size(cl));
h5write(savename, '/RAA/CL', cl);

h5create(savename,'/RAA/mN41Sq_a',size(f_RAAa(:,2)));
h5write(savename, '/RAA/mN41Sq_a', f_RAAa(:,2))
h5create(savename,'/RAA/sin2T4Sq_a',size(f_RAAa(:,2)));
h5write(savename, '/RAA/sin2T4Sq_a', f_RAAa(:,1))

h5create(savename,'/RAA/mN41Sq_b',size(f_RAAb(:,2)));
h5write(savename, '/RAA/mN41Sq_b', f_RAAb(:,2))
h5create(savename,'/RAA/sin2T4Sq_b',size(f_RAAb(:,2)));
h5write(savename, '/RAA/sin2T4Sq_b', f_RAAb(:,1))

h5create(savename,'/RAA/mN41Sq_c',size(f_RAAc(:,2)));
h5write(savename, '/RAA/mN41Sq_c', f_RAAc(:,2))
h5create(savename,'/RAA/sin2T4Sq_c',size(f_RAAc(:,2)));
h5write(savename, '/RAA/sin2T4Sq_c', f_RAAc(:,1))

%% BEST and Gallium

f_BestGa1a = importdata([fdir,'contour_BestGAcombi_1sigma_a.txt']);
f_BestGa1b = importdata([fdir,'contour_BestGAcombi_1sigma_b.txt']);
f_BestGa2a = importdata([fdir,'contour_BestGAcombi_2sigma_a.txt']);
f_BestGa2b = importdata([fdir,'contour_BestGAcombi_2sigma_b.txt']);
f_BestGa3 = importdata([fdir,'contour_BestGAcombi_3sigma.txt']);

str = "https://arxiv.org/abs/2109.11482";
h5create(savename,'/BestGaCombi/Reference',size(str),'Datatype','string');
h5write(savename,'/BestGaCombi/Reference',str)
h5create(savename,'/BestGaCombi/1Sigma/mN41Sq_a',[numel(f_BestGa1a(:,2)),1]);
h5write(savename, '/BestGaCombi/1Sigma/mN41Sq_a', f_BestGa1a(:,2));
h5create(savename,'/BestGaCombi/1Sigma/mN41Sq_b',[numel(f_BestGa1b(:,2)),1]);
h5write(savename, '/BestGaCombi/1Sigma/mN41Sq_b', f_BestGa1b(:,2));
h5create(savename,'/BestGaCombi/1Sigma/sin2T4Sq_a',[numel(f_BestGa1a(:,1)),1]);
h5write(savename, '/BestGaCombi/1Sigma/sin2T4Sq_a', f_BestGa1a(:,1));
h5create(savename,'/BestGaCombi/1Sigma/sin2T4Sq_b',[numel(f_BestGa1b(:,1)),1]);
h5write(savename, '/BestGaCombi/1Sigma/sin2T4Sq_b', f_BestGa1b(:,1));

h5create(savename,'/BestGaCombi/2Sigma/mN41Sq_a',[numel(f_BestGa2a(:,2)),1]);
h5write(savename, '/BestGaCombi/2Sigma/mN41Sq_a', f_BestGa2a(:,2));
h5create(savename,'/BestGaCombi/2Sigma/mN41Sq_b',[numel(f_BestGa2b(:,2)),1]);
h5write(savename, '/BestGaCombi/2Sigma/mN41Sq_b', f_BestGa2b(:,2));
h5create(savename,'/BestGaCombi/2Sigma/sin2T4Sq_a',[numel(f_BestGa2a(:,1)),1]);
h5write(savename, '/BestGaCombi/2Sigma/sin2T4Sq_a', f_BestGa2a(:,1));
h5create(savename,'/BestGaCombi/2Sigma/sin2T4Sq_b',[numel(f_BestGa2b(:,1)),1]);
h5write(savename, '/BestGaCombi/2Sigma/sin2T4Sq_b', f_BestGa2b(:,1));

h5create(savename,'/BestGaCombi/3Sigma/mN41Sq',[numel(f_BestGa3(:,2)),1]);
h5write(savename, '/BestGaCombi/3Sigma/mN41Sq', f_BestGa3(:,2));
h5create(savename,'/BestGaCombi/3Sigma/sin2T4Sq',[numel(f_BestGa3(:,1)),1]);
h5write(savename, '/BestGaCombi/3Sigma/sin2T4Sq', f_BestGa3(:,1));

%% DANSS
f_DANSS = importdata([fdir,'coord_DANSS_95CL.mat']);
h5create(savename,'/DANSS/Reference',size(string(f_DANSS.Reference)),'Datatype','string');
h5write(savename,'/DANSS/Reference',string(f_DANSS.Reference))
h5create(savename,'/DANSS/mN41Sq',size(f_DANSS.DmSquare41_Y));
h5write(savename, '/DANSS/mN41Sq', f_DANSS.DmSquare41_Y)
h5create(savename,'/DANSS/sin2T4Sq',size(f_DANSS.SinSquare2Theta_X));
h5write(savename, '/DANSS/sin2T4Sq', f_DANSS.SinSquare2Theta_X)
h5create(savename,'/DANSS/CL',size(str2double(f_DANSS.CL)));
h5write(savename, '/DANSS/CL', str2double(f_DANSS.CL));

%% PROSPECT
f_PROSPECT = importdata([fdir,'coord_Prospect2020_95CL.mat']);
h5create(savename,'/Prospect/Reference',size(string(f_PROSPECT.Reference)),'Datatype','string');
h5write(savename,'/Prospect/Reference',string(f_PROSPECT.Reference))
h5create(savename,'/Prospect/mN41Sq',size(f_PROSPECT.DmSquare41_Y));
h5write(savename, '/Prospect/mN41Sq', f_PROSPECT.DmSquare41_Y)
h5create(savename,'/Prospect/sin2T4Sq',size(f_PROSPECT.SinSquare2Theta_X));
h5write(savename, '/Prospect/sin2T4Sq', f_PROSPECT.SinSquare2Theta_X)
h5create(savename,'/Prospect/CL',size(str2double(f_PROSPECT.CL)));
h5write(savename, '/Prospect/CL', str2double(f_PROSPECT.CL));

%% STEREO
f_Stereo = importdata([fdir,'coord_STEREOprd102_95CL.mat']);
h5create(savename,'/Stereo/Reference',size(string(f_Stereo.Reference)),'Datatype','string');
h5write(savename,'/Stereo/Reference',string(f_Stereo.Reference))
h5create(savename,'/Stereo/mN41Sq',size(f_Stereo.DmSquare41_Y));
h5write(savename, '/Stereo/mN41Sq', f_Stereo.DmSquare41_Y)
h5create(savename,'/Stereo/sin2T4Sq',size(f_Stereo.SinSquare2Theta_X));
h5write(savename, '/Stereo/sin2T4Sq', f_Stereo.SinSquare2Theta_X)
h5create(savename,'/Stereo/CL',size(str2double(f_Stereo.CL)));
h5write(savename, '/Stereo/CL', str2double(f_Stereo.CL));

%% Neutrino 4
f_N4 = importdata([fdir,'coord_Neutrino4_123sigma.mat']);
h5create(savename,'/Neutrino4/Reference',size(string(f_N4.Reference)),'Datatype','string');
h5write(savename,'/Neutrino4/Reference',string(f_N4.Reference));
h5create(savename,'/Neutrino4/1Sigma/mN41Sq',size(f_N4.DmSquare41_Y_1sigma));
h5write(savename, '/Neutrino4/1Sigma/mN41Sq', f_N4.DmSquare41_Y_1sigma);
h5create(savename,'/Neutrino4/1Sigma/sin2T4Sq',size(f_N4.SinSquare2Theta_X_1sigma));
h5write(savename, '/Neutrino4/1Sigma/sin2T4Sq', f_N4.SinSquare2Theta_X_1sigma);

h5create(savename,'/Neutrino4/2Sigma/mN41Sq',size(f_N4.DmSquare41_Y_2sigma));
h5write(savename, '/Neutrino4/2Sigma/mN41Sq', f_N4.DmSquare41_Y_2sigma);
h5create(savename,'/Neutrino4/2Sigma/sin2T4Sq',size(f_N4.SinSquare2Theta_X_2sigma));
h5write(savename, '/Neutrino4/2Sigma/sin2T4Sq', f_N4.SinSquare2Theta_X_2sigma);

h5create(savename,'/Neutrino4/3Sigma/mN41Sq',size(f_N4.DmSquare41_Y_3sigma));
h5write(savename, '/Neutrino4/3Sigma/mN41Sq', f_N4.DmSquare41_Y_3sigma);
h5create(savename,'/Neutrino4/3Sigma/sin2T4Sq',size(f_N4.SinSquare2Theta_X_3sigma));
h5write(savename, '/Neutrino4/3Sigma/sin2T4Sq', f_N4.SinSquare2Theta_X_3sigma);

%% Mainz
f_Mainz = importdata([fdir,'coord_Mainz_95CL.mat']);
h5create(savename,'/Mainz/Reference',size(string(f_Mainz.Reference)),'Datatype','string');
h5write(savename,'/Mainz/Reference',string(f_Mainz.Reference))
h5create(savename,'/Mainz/mN41Sq',size(f_Mainz.DmSquare41_Y));
h5write(savename, '/Mainz/mN41Sq', f_Mainz.DmSquare41_Y)
h5create(savename,'/Mainz/sin2T4Sq',size(f_Mainz.SinSquare2Theta_X));
h5write(savename, '/Mainz/sin2T4Sq', f_Mainz.SinSquare2Theta_X)
h5create(savename,'/Mainz/CL',size(str2double(f_Mainz.CL)));
h5write(savename, '/Mainz/CL', str2double(f_Mainz.CL));

%% Troitsk
f_Troitsk = importdata([fdir,'coord_Troitsk_95CL.mat']);
h5create(savename,'/Troitsk/Reference',size(string(f_Troitsk.Reference)),'Datatype','string');
h5write(savename,'/Troitsk/Reference',string(f_Troitsk.Reference))
h5create(savename,'/Troitsk/mN41Sq',size(f_Troitsk.DmSquare41_Y));
h5write(savename, '/Troitsk/mN41Sq', f_Troitsk.DmSquare41_Y)
h5create(savename,'/Troitsk/sin2T4Sq',size(f_Troitsk.SinSquare2Theta_X));
h5write(savename, '/Troitsk/sin2T4Sq', f_Troitsk.SinSquare2Theta_X)
h5create(savename,'/Troitsk/CL',size(str2double(f_Troitsk.CL)));
h5write(savename, '/Troitsk/CL', str2double(f_Troitsk.CL));




