function target_distribution = f_s2s_wave_propergation_GPU( source_distribution, ...
    focus_length, wave_length, lens2source, target2lens, lens_radius, ...
    target_distribution_init)
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
% check to ensure the values are in range
% assert(mod(sz_source_distribution(1),2)==0,mod(sz_source_distribution(2),2)==0,'Each side of the source plain should not contain odd number of pixels!')
assert(focus_length>0,'focus_length should not be equal nor smaller than 0!')
assert(wave_length>0,'wave_length should not be equal nor smaller than 0!')
assert(lens2source>0,'lens2source should not be equal nor smaller than 0!')
assert(target2lens>0,'target2lens should not be equal nor smaller than 0!')
assert(lens_radius>0,'lens_radius should not be equal nor smaller than 0!')
assert(lens_radius>0,'target should not be at the same side with the source!')
%% Calculations
%% Initialize
target_distribution = target_distribution_init;
%% simple notations
Us = source_distribution{1}; % all values of the source plain
res_s = source_distribution{2}; % The resolution of the pixel
h_s = sz_source_distribution(1); % The hight of the source plain in pixels
w_s = sz_source_distribution(2); % The width of the source plain in pixels

f = focus_length;
l = wave_length;
a = lens2source;
b = target2lens;
r = lens_radius;
zt = lens2source + target2lens;

k = 2*pi/l;
alpha = -pi/2*(1/a+1/b-1/f);
assert(abs(alpha)>1e-6,'target plain should not be near the image plain!') 

%% preparations
% The limits for the integration field
y_s_min = -h_s/2*res_s;
y_s_max = h_s/2*res_s;
x_s_min = -w_s/2*res_s;
x_s_max = w_s/2*res_s;
disp(x_s_min)
disp(x_s_max)
disp(y_s_min)
disp(y_s_max)
% Build meshgrid
x_s = x_s_min:(x_s_max-x_s_min)/(w_s-1):x_s_max;
y_s = y_s_min:(y_s_max-y_s_min)/(h_s-1):y_s_max;
[Xs,Ys] = meshgrid(x_s,y_s);
Xs_GPU = gpuArray(Xs);
Ys_GPU = gpuArray(Ys);
Us_GPU = gpuArray(Us);
% Calculate the phase angle if the electical field: E = u(x,y)*exp(1i*theta)
f_theta = @(yt,xt)-k/2*((Xs_GPU.^2+Ys_GPU.^2)/a+(xt^2+yt^2)/b+k/2/alpha*(Xs_GPU/a+xt/b).^2+k/2/alpha*(Ys_GPU/a+yt/b).^2);
% Calculate the integration kernal: E = u(x,y)*exp(1i*theta)
f_kernal = @(yt,xt)Us_GPU.*exp(1i*f_theta(yt,xt));
% A constant scaling factor for the integration
const = pi*exp(-1i*3/2*pi)/l^2/a/b*exp(-1i*k*zt)*exp(-1i*k/2/f*r^2)/alpha;

%% Calculations
td_size = size(target_distribution{1});
td_h = td_size(1);
td_w = td_size(2);
td_res = target_distribution{2};
td = gpuArray(target_distribution{1});
% % The limits for the target field
% y_t_min = -td_h/2*td_res;
% y_t_max = td_h/2*td_res;
% x_t_min = -td_w/2*td_res;
% x_t_max = td_w/2*td_res;
% % Build meshgrid
% x_t = x_t_min:(x_t_max-x_t_min)/(td_w-1):x_t_max;
% y_t = y_t_min:(y_t_max-y_t_min)/(td_h-1):y_t_max;
% [Xt,Yt] = meshgrid(x_t,y_t);
% Xt_GPU = gpuArray(Xt);
% Yt_GPU = gpuArray(Yt);

y_s_GPU = gpuArray(y_s);
x_s_GPU = gpuArray(x_s);

% td = arrayfun(@(yt,xt)const*trapz(y_s_GPU,trapz(x_s_GPU,f_kernal(yt,xt),2)),Yt_GPU,Xt_GPU);

for row = 1:td_h
    for col = 1:td_w
        xt = (col-td_w/2)*td_res;
        yt = (row-td_h/2)*td_res;
        % integration
        F = f_kernal(yt,xt);
        td(row,col) = const*trapz(y_s_GPU,trapz(x_s_GPU,F,2));
    end
end
target_distribution{1} = gather(td);
end



