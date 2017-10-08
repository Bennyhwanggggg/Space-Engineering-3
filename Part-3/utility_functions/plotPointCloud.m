function plotPointCloud(fNum,fName)
    f = figure(fNum);
    f.Name = fName;
    load('localTerrainPointCloud.mat');
    x = pointCloud(:,1);
    y = pointCloud(:,2);
    z = pointCloud(:,3);
    plot3(x,y,z,'k.');
    xlabel('[m]');
    ylabel('[m]');
    zlabel('[m]');
    axis equal
    view(170,24);
end