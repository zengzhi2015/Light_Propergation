close all
clear
clc
%% Set parameters
focus_length = 2e-2;
wave_length = 1e-14;
lens2source = 5e-2;
target2lens = 3e-2; 
lens_radius = 1e-2;

source_distribution = {zeros(512,512),2e-2/512};
target_position = [0,0];
%% Checkings
% check the shape of the inputs
sz_source_distribution = size(source_distribution);
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
% simple notations
Us = source_distribution{1};
res = source_distribution{2};
h = sz_source_distribution(1);
w = sz_source_distribution(2);

f = focus_length;
l = wave_length;
a = lens2source;
b = target2lens;
r = lens_radius;
xt = target_position(1);
yt = target_position(2);
zt = lens2source + target2lens;

%short note
k = 2*pi/l;
alpha = -pi/2*(1/a+1/b-1/f);
assert(abs(alpha)>1e-6,'target plain should not be near the image plain!') 

%preparations
f_us = @(ys,xs)Us(max(1,min(h,round(h/2+ys/res))),max(1,min(h,round(w/2+xs/res))));
f_theta = @(ys,xs)-k/2*((xs.^2+ys.^2)/a+(xt^2+yt^2)/b+k/2/alpha*(xs/a+xt/b).^2+k/2/alpha*(ys/a+yt/b).^2);
f_kernal = @(ys,xs)f_us(ys,xs).*exp(1i*f_theta(ys,xs));
const = pi*exp(-1i*3/2*pi)/l^2/a/b*exp(-1i*k*zt)*exp(-1i*k/2/f*r^2)/alpha;
ymin = -h/2*res;
ymax = h/2*res;
xmin = -w/2*res;
xmax = w/2*res;