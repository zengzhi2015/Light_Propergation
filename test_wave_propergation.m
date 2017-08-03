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
%% Plot the lens
THETA = 0:0.01:2*pi;
X = lens_radius*cos(THETA);
Y = lens_radius*sin(THETA);
Z = 0*X + lens2source;
C = X*0 + 0.5;
fill3(X,Y,Z,C)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
% Plot the source plane
hold on
sp_sz = size(source_distribution{1});
res = source_distribution{2};
X = [sp_sz(2)/2*res, sp_sz(2)/2*res, -sp_sz(2)/2*res, -sp_sz(2)/2*res, sp_sz(2)/2*res];
Y = [sp_sz(1)/2*res, -sp_sz(1)/2*res, -sp_sz(1)/2*res, sp_sz(1)/2*res, sp_sz(1)/2*res];
Z = 0*X;
C = X*0 + 0.3;
fill3(X,Y,Z,C)
% Plot the target plane
hold on
sp_sz = size(source_distribution{1});
res = source_distribution{2};
X = [sp_sz(2)/2*res, sp_sz(2)/2*res, -sp_sz(2)/2*res, -sp_sz(2)/2*res, sp_sz(2)/2*res];
Y = [sp_sz(1)/2*res, -sp_sz(1)/2*res, -sp_sz(1)/2*res, sp_sz(1)/2*res, sp_sz(1)/2*res];
Z = 0*X + lens2source + target2lens;
C = X*0 + 0.7;
fill3(X,Y,Z,C)
%% Calculations
target_distribution = f_wave_propergation( source_distribution, ...
    focus_length, wave_length, lens2source, target2lens, lens_radius, target_position);