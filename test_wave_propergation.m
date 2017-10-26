close all
clear
clc
%% Set parameters
focus_length = 50e-3;
wave_length = 2e-3;
lens2source = 50e-3;
target2lens = 30e-3; 
lens_radius = 10e-3;
% The following parameters are for debugging
% focus_length = 2;
% wave_length = 0.02;
% lens2source = 5;
% target2lens = 3; 
% lens_radius = 1;

source_distribution = {zeros(100,100)+0.01,10e-3/1024};
target_distribution = {zeros(100,100),10e-3/1024};
% source_distribution = {zeros(400,400)+1,0.64/512};
% target_distribution = {zeros(100,100),6.4/512};
% target_position = [0,0];

%% Calculations
td_size = size(target_distribution{1});
td_h = td_size(1);
td_w = td_size(2);
td_res = target_distribution{2};
td = target_distribution{1};

parfor (row = 1:td_h, 8)
%parfor row = 1:td_h
    for col = 1:td_w
        td(row,col) = f_s2p_wave_propergation( source_distribution, ...
            focus_length, wave_length, lens2source, target2lens, lens_radius, ...
            [(col-td_w/2)*td_res,(row-td_h/2)*td_res]);
    end
end
target_distribution{1} = td;
% target_value = f_s2p_wave_propergation( source_distribution, ...
%     focus_length, wave_length, lens2source, target2lens, lens_radius, target_position);
%% Plot distributions
figure(1)
surf(abs(td),'EdgeAlpha',0.2)
%% Plot the lens
figure(2)
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
sp_sz = size(target_distribution{1});
res = target_distribution{2};
X = [sp_sz(2)/2*res, sp_sz(2)/2*res, -sp_sz(2)/2*res, -sp_sz(2)/2*res, sp_sz(2)/2*res];
Y = [sp_sz(1)/2*res, -sp_sz(1)/2*res, -sp_sz(1)/2*res, sp_sz(1)/2*res, sp_sz(1)/2*res];
Z = 0*X + lens2source + target2lens;
C = X*0 + 0.7;
fill3(X,Y,Z,C)