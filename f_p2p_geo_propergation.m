function [target, interm] = f_p2p_geo_propergation( focus_length, point_source, target_z )
%F_P2P_PROPERGATION Summary of this function goes here
% This function calculate the projection of a point light source to the
% target plain whose z value is target_z
%   Detailed explanation goes here
%   The coordinate id defined as follows:
%       O: at the center of the lens
%       Z: from the center of the lens toward the space of sources and
%       perpendicular to the lens plain
%       X: in the lens plain and perpendicular to the Z axis
%       Y: perpendicular to X and Z, and follow the right hand rule
%   Inputs:
%       point_source: [xs,ys,zs,dxs,dys,dzs,intensity]
%       target_z: zt
%   Outputs:
%       target: [xt,yt,xt,dxt,dyt,dzt,intensity]

%% Checkings

% check the shape of the inputs
sz_focus_length = size(focus_length);
assert(sz_focus_length(1)==1 && sz_focus_length(2)==1 && length(sz_focus_length)==2,'The shape of focus_length is not correct!')
sz_point_source = size(point_source);
assert(sz_point_source(1)==1 && sz_point_source(2)==7 && length(sz_point_source)==2,'The shape of sz_point_source is not correct!')
sz_target_z = size(target_z);
assert(sz_target_z(1)==1 && sz_target_z(2)==1 && length(sz_target_z)==2,'The shape of sz_target_z is not correct!')

% check to ensure that focus_length is greater than zero
assert(focus_length>0,'focus_length should not be equal nor smaller than 0!')

xs = point_source(1);
ys = point_source(2);
zs = point_source(3);
dxs = point_source(4);
dys = point_source(5);
dzs = point_source(6);
intensity = point_source(7);

% check that the z position of point source is greater than zero
assert(zs>0,'the z position of point source should not be equal nor smaller than 0!')

% check that the dz component of the direction of point source is smaller than zero
assert(dzs<0,'the dz position of point source should not be equal nor smaller than 0!')

% check to ensure that target_z is smaller than zero
assert(target_z<0,'target_z should not be equal nor larger than 0!')

%% Calculations

% Calculate the distance from the point source to the lens plain:
source2lens = abs(zs);

% calculate the projected point on the lens plain
xl = xs - source2lens / dzs * dxs;
yl = ys - source2lens / dzs * dys;
zl = 0;

% assert(-(xl-xs)/source2lens == dxs/dzs && -(yl-ys)/source2lens == dys/dzs && -(zl-zs)/source2lens == 1, ... 
%     'The calculation of the projection to lens plain is not correct!')

% calculate the projected point on the focal plain where a ray comes from
% the origin and toward the direction which is the same as that of the
% point source. According to the theory of parrllel lightening, this is the
% point that all light that parrallel to that of the source would focus on.

xf = 0 - focus_length / dzs * dxs;
yf = 0 - focus_length / dzs * dys;
zf = -focus_length;

assert(-xf/focus_length == dxs/dzs && -yf/focus_length == dys/dzs && -zf/focus_length == 1, ... 
    'The calculation of the projection to lenze plain is not correct!')

% Calculate the new direction on the lens, focal, or the target plain
dxt = xf - xl;
dyt = yf - yl;
dzt = zf - zl;

assert(norm([dxt,dyt,dzt])>0, ...
    'the distance from the lens plain point to the focal plain point is wrong!')

norm_temp = norm([dxt,dyt,dzt]);
dxt = dxt / norm_temp;
dyt = dyt / norm_temp;
dzt = dzt / norm_temp;

% assert(norm([dxt,dyt,dzt])==1, ...
%     'the direction of target point is not correct!')

% Calculate the target point
target2lens = abs(target_z);
xt = xl - dxt / dzt * target2lens;
yt = yl - dyt / dzt * target2lens;
zt = zl - target2lens;

% assert(-(xt-xl)/target2lens == dxt/dzt && -(yt-yl)/target2lens == dyt/dzt && -(zt-zl)/target2lens == 1, ... 
%     'The calculation of the target is not correct!')

target = [xt,yt,zt,dxt,dyt,dzt,intensity];
interm = [xl,yl,zl];

%% Plot the result for checking
% hold on
% plot3([xs,xl],[ys,yl],[zs,zl],'g')
% plot3([xl,xt],[yl,yt],[zl,zt],'r')
% 
% scatter3([xs,xt],[ys,yt],[zs,zt],[10,10],[intensity,intensity],'filled')

end

