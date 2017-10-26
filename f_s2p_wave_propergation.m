function target_distribution = f_s2p_wave_propergation( source_distribution, ...
    focus_length, wave_length, lens2source, target2lens, lens_radius, target_position)
%F_WAVE_PROPERGATION Summary of this function goes here
% This function is used to calculate the electric field propergation of
% light using wave oprics theory. Given the distribution of the
% electric field on a plain parrallel to the lens, we calculate the
% electric field of an arbitraty point in the opposite side of the lens,
% except for those on the image plain (1/a + 1/b != 1/f).
%   Detailed explanation goes here
%   The coordinate id defined as follows:
%       O: at the center of source plain
%       Z: from the center of the source plain toward that of the lens
%       X: in the lens plain and perpendicular to the Z axis
%       Y: perpendicular to X and Z, and follow the right hand rule
%   Inputs:
%       source_distribution: {[us, ...;us, ...; ...],resolution}
%       focus_length: f
%       wave_length: l
%       lens2source: a
%       target2lens: b
%       lens_radius: r
%       target_position: [xt,yt]
%   Outputs:
%       target_distribution: ut
%   Note:
%       

%% Checkings
% check the shape of the inputs
sz_source_distribution = size(source_distribution{1});
assert(length(sz_source_distribution)==2,'The shape of source_distribution is not correct!')
sz_focus_length = size(focus_length);
assert(sz_focus_length(1)==1 && sz_focus_length(2)==1 && length(sz_focus_length)==2,'The shape of focus_length is not correct!')
sz_wave_length = size(wave_length);
assert(sz_wave_length(1)==1 && sz_wave_length(2)==1 && length(sz_wave_length)==2,'The shape of wave_length is not correct!')
sz_lens2source = size(lens2source);
assert(sz_lens2source(1)==1 && sz_lens2source(2)==1 && length(sz_lens2source)==2,'The shape of lens2source is not correct!')
sz_target2lens = size(target2lens);
assert(sz_target2lens(1)==1 && sz_target2lens(2)==1 && length(sz_target2lens)==2,'The shape of target2lens is not correct!')
sz_lens_radius = size(lens_radius);
assert(sz_lens_radius(1)==1 && sz_lens_radius(2)==1 && length(sz_lens_radius)==2,'The shape of lens_radius is not correct!')
sz_target_position = size(target_position);
assert(sz_target_position(1)==1 && sz_target_position(2)==2 && length(sz_focus_length)==2,'The shape of target_position is not correct!')
% check to ensure the values are in range
% assert(mod(sz_source_distribution(1),2)==0,mod(sz_source_distribution(2),2)==0,'Each side of the source plain should not contain odd number of pixels!')
assert(focus_length>0,'focus_length should not be equal nor smaller than 0!')
assert(wave_length>0,'wave_length should not be equal nor smaller than 0!')
assert(lens2source>0,'lens2source should not be equal nor smaller than 0!')
assert(target2lens>0,'target2lens should not be equal nor smaller than 0!')
assert(lens_radius>0,'lens_radius should not be equal nor smaller than 0!')
assert(lens_radius>0,'target should not be at the same side with the source!')
%% Calculations
%% simple notations
Us = source_distribution{1}; % all positions of the source plain
res = source_distribution{2}; % The resolution of the pixel
h = sz_source_distribution(1); % The hight of the source plain in pixels
w = sz_source_distribution(2); % The width of the source plain in pixels

f = focus_length;
l = wave_length;
a = lens2source;
b = target2lens;
r = lens_radius;
xt = target_position(1);
yt = target_position(2);
zt = lens2source + target2lens;

%% short note
k = 2*pi/l;
alpha = -pi/2*(1/a+1/b-1/f);
assert(abs(alpha)>1e-6,'target plain should not be near the image plain!') 

%% preparations
% (m,m) -> (pixels,pixels) -> electral field intensity at that position
f_us = @(ys,xs)Us(max(1,min(h,round(h/2+ys/res)))+(max(1,min(h,round(w/2+xs/res)))-1)*h);
% Calculate the phase angle if the electical field: E = u(x,y)*exp(1i*theta)
f_theta = @(ys,xs)-k/2*((xs.^2+ys.^2)/a+(xt^2+yt^2)/b+k/2/alpha*(xs/a+xt/b).^2+k/2/alpha*(ys/a+yt/b).^2);
% Calculate the integration kernal: E = u(x,y)*exp(1i*theta)
f_kernal = @(ys,xs)f_us(ys,xs).*exp(1i*f_theta(ys,xs));
% A constant scaling factor for the integration
const = pi*exp(-1i*3/2*pi)/l^2/a/b*exp(-1i*k*zt)*exp(-1i*k/2/f*r^2)/alpha;
% The limits for the integration field
ymin = -h/2*res;
ymax = h/2*res;
xmin = -w/2*res;
xmax = w/2*res;

% integration
% target_distribution = const*integral2(f_kernal,ymin,ymax,xmin,xmax,'RelTol',1e-12,'AbsTol',1e-12);
% target_distribution = const*integral2(f_kernal,ymin,ymax,xmin,xmax);
% target_distribution = const*quad2d(f_kernal,ymin,ymax,xmin,xmax);
x = xmin:(xmax-xmin)/100:xmax;
y = ymin:(ymax-ymin)/100:ymax;
[X,Y] = meshgrid(x,y);
F = f_kernal(X,Y);
target_distribution = const*trapz(y,trapz(x,F,2));
end

