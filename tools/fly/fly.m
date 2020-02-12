function varargout= fly(varargin)
%FLY fly through a 3D figure
%	FLY allows to freely move through any 3D figure piloting with the keyboard.
%	The figure has to be generated before activating the function.
%	Changes in direction or speed are maintained until compensated by the
%	opposite keystroke. Repeated strokes of a key accelerate the
%	corresponding movement.
%
%	To pilot use:
%		Right, Left, Up or Down: arrow keys
%		Accelerate forward:		'a'
%		Brake or move backward: 'z'
%		Full stop: backspace
%		Quit: escape
%	Controls can be changed in function 'key_stroke'
% 
%	fly('start',V) Positions the camera and target at the start at coordinates 
%	specified in vector V. V must be a 6-element vector. The first three values
%	are the X, Y, Z coordinates of the camera. The last three values are
%	the coordinates of the target (the direction you are looking to).
%	If not specified a starting point is chosen by default. 
%
%	flight= FLY records and returns the coordinates of camera and target 
%	during the flight as Nx6 array 'flight'. N is the number of frames. The first
%	three columns are the x,y,z coordinates of the camera, and the last three 
%	those of the target. This array can be used to generate a video of the 
%	flight with another function (see template function flyvideo.m).
%	No recording is made while you are not moving.
%
%	An example figure is provided in fly_figure.m taken from MATLAB
%	documentation. A video example is also included in the submission.
%
% 	Example:
% 	fly_figure
% 	fly
%
% Author: Francisco de Castro (2014)

global elev azim speed keepflying step


% Defaults
record= 'on';
start= [];


% Optional parameters
for j= 1:2:size(varargin,2)
	switch varargin{j}
		case 'start'
			start= varargin{j+1};
			if numel(start) < 6
				error('Starting location vector must have at least 6 elements'); 
			end
		otherwise
			error('Unknown parameter. Accepted values: show, plotfit')
	end
end


% Save flight?
if nargout > 0, record= 'on'; end


% Get current figure
figure(gcf)


% Must be 3D stuff
view(3)
axis vis3d


% Keep flying by default
keepflying= 1;


% Get keystrokes
set (gcf, 'KeyPressFcn', @key_stroke);


% Get coordinate limits and size
xa= xlim; ya= ylim; za= zlim;
diagonal= sqrt(range(xa)^2 + range(ya)^2 + range(za)^2);


% If no starting poin, start from one corner looking to the oppostite
if isempty(start)
	C= [xa(2) ya(2) za(2)];
	T= [xa(1) ya(1) za(1)];
else
	C= start(1:3);
	T= start(4:6);
end


% Initial speed of camera and target (both static) and acceleration (linear
% and angular)
elev=  0;
azim=  0;
speed= 0;
step= [diagonal*1e-5 5e-5];


% Keep flying until 'esc' is pressed
flypath= [C(1:3), T(1:3)];
while keepflying
	[C,T]= move_position(T,C,azim,elev,speed,diagonal);
	campos(C);
	camtarget(T);
	drawnow;
	if strcmp(record,'on') && speed ~= 0 && elev ~= 0 && azim ~= 0
		flypath= [flypath; C(1:3), T(1:3)];
	end
end


% Out of fuel
close(gcf)
if nargout > 0 varargout= {flypath}; end
end



function [C,T]= move_position(T,C,azim,elev,speed,diagonal)

	% Rotate target around camera
	[th,phi,r]= cart2sph(T(1)-C(1), T(2)-C(2), T(3)-C(3));
	[X(1), X(2), X(3)]= sph2cart(th + azim, phi + elev, r);
	T= X + C;

%	Keep target inside box? Maybe necessary to avoid figure from dissappearing
% 	T= [max(min(T(1),xa(2)),xa(1)), max(min(T(2),ya(2)),ya(1)), max(min(T(3),za(2)),za(1))];

	% Close-in camera to target
	[th,phi,r]= cart2sph(C(1)-T(1), C(2)-T(2), C(3)-T(3));
	dist= r-speed;
	[X(1), X(2), X(3)]= sph2cart(th, phi, dist);
	C= X + T;

	% If too close, push target back a little
	if dist < diagonal/500
		[th,phi,r]= cart2sph(T(1)-C(1), T(2)-C(2), T(3)-C(3));
		[X(1), X(2), X(3)]= sph2cart(th, phi, diagonal/500);
		T= X + C;
	end
end



function key_stroke(src,eventdata)
global elev azim speed keepflying step

% Move depending on key stroke
switch eventdata.Key

	case 'a' % Accelerate forward
		speed= speed+ step(1);
	case 'z' % Brake and backwards
		speed= speed- step(1);

	case 'leftarrow' % Turn left
		azim= azim+ step(2);

	case 'rightarrow' % Turn right
		azim= azim- step(2);

	case 'uparrow' % Look down (pitch)
		elev= elev- step(2);

	case 'downarrow' % Look up (pitch)
		elev= elev+ step(2);

	% Full stop
	case 'backspace'
		speed= 0; elev= 0; azim= 0;

	% Tired of flying
	case 'escape'
		keepflying= 0;

end

end
