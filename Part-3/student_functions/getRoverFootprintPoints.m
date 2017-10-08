%% function find rover foot print points
%
% Author Kuan Chun Hwang

function points = getRoverFootprintPoints(pos,ROVER_FOOTPRINT_RADIUS)

% Find cloud points and store in pointCloudps
load('localTerrainPointCloud.mat');
x = pointCloud(:,1);
y = pointCloud(:,2);
z = pointCloud(:,3);
pointCloudpos = [x y];

% Search for indext number for the points that are inside the radius
idx = rangesearch(pointCloudpos,pos,ROVER_FOOTPRINT_RADIUS);

% Find the points
points = pointCloud(idx{1,:},:);
end