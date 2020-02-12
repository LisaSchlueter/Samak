function flyvideo(filename,flypath,framestep)
% FLYVIDEO Creates a video from output of fly.m
%	This is intended as a template function that can be modified 
%  to accomodate specific needs.
% 
% filename: name of the file to store the video. AVI format by default.
% flypath: array of size Nx6. N is the number of steps in the flight. The
%		first 3 columns are the x, y, z coordinates of the camera and the
%		last 3 columns are the same for the target.
% framestep: one out of framestep frames is saved into the video. Lower
%		numbers mean higher quality and larger files.
% 
% Before using this function you must recreate the *same* figure used to generate
% the flight.
% 
%	Author: Francisco de Castro

% Open video file
vidObj = VideoWriter(filename);
open(vidObj);

% Get figure and set perspective
figure(gcf)
view(3)
axis vis3d
set(gcf,'renderer','zbuffer'); % Otherwise writeVideo doesn't work

% Reproduce flight. 
for j= 1:framestep:size(flypath,1)
	campos(flypath(j,1:3));
	camtarget(flypath(j,4:6));
	drawnow;
	currFrame = getframe(gcf);
	writeVideo(vidObj,currFrame);
end

% Close all
close(vidObj)
close(gcf)
