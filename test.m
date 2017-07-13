close all
clear
clc
%%
THETA = 0:0.01:2*pi;
X = 2*cos(THETA);
Y = 2*sin(THETA);
Z = 0*X;
C = X*0 + 0.5;
fill3(X,Y,Z,C)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
%%
focus_length = 1;

% generate 1000 random points
number_test = 1000;
point_sources = zeros([number_test,7]);
targets = zeros([number_test,7]);
interms = zeros([number_test,3]);

for i = 1:2000
    r = 0.5*rand();
    theta = 2*pi*rand();
    xs = r*cos(theta);
    ys = r*sin(theta);
    zs = 4;
    dxs = 0.2 - xs*0.05;
    dys = 0.1 + ys*0.15;
    dzs = -1;
    intensity = exp(-r^2);
    point_sources(i,:) = [xs,ys,zs,dxs,dys,dzs,intensity];
    target_z = -2;
    [targets(i,:),interms(i,:)] = f_p2p_geo_propergation( focus_length, point_sources(i,:), target_z );
end

%%
hold on

plot3([point_sources(:,1),interms(:,1)]', ...
    [point_sources(:,2),interms(:,2)]', ...
    [point_sources(:,3),interms(:,3)]','g')

plot3([interms(:,1),targets(:,1)]', ...
    [interms(:,2),targets(:,2)]', ...
    [interms(:,3),targets(:,3)]','r')

X2 = [point_sources(:,1)',targets(:,1)'];
Y2 = [point_sources(:,2)',targets(:,2)'];
Z2 = [point_sources(:,3)',targets(:,3)'];
S2 = X2*0 + 10;
C2 = [point_sources(:,7)',targets(:,7)'];
scatter3(X2,Y2,Z2,S2,C2,'filled')



