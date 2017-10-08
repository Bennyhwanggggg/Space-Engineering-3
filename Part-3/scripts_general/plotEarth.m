%%  Plot Earth
%
% Author Kuan Chun Hwang

% load topographical Earth map
load('topo.mat','topo');

% Create a sphere, make it earth sized (in meters)
[x,y,z] = sphere(50);
x = -x.*6378000;
y = -y.*6378000;
z = z.*6378000;

props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;

% Plot Earth
axes('dataaspectratio',[1 1 1],'visible','on')
title('3D ECI plot');
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
hold on
sh=surface(x,y,z,props);