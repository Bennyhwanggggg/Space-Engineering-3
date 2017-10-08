%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlotEarthSphere.m - Script for plotting an Earth sized sphere to use for
% visualisation of you satellite orbit. You
% may either use this script or cut and paste the code into your own
% script.
%
% By Kuan-Chun Hwang.
%
close all
clear
clc

% Create figure and load topographical Earth map
figure
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
hold on
sh=surface(x,y,z,props);

% Time stamp per step (s)
t_stamp = 60*60;

for i=1:864000
    rotate(sh,[0,0,1],(t_stamp*7.292115e-5));
    pause(0.01);
end

