function [x,y,z] = sampleMap(numPoints)
    
    load('localTerrainRed.mat');
    
    gridSize = length(xm);
    
    xIndex = randi(gridSize,numPoints);
    yIndex = randi(gridSize,numPoints);
    x = zeros(1,numPoints);
    y = zeros(1,numPoints);
    z = zeros(1,numPoints);
    
    for i = 1:numPoints
        z(i) = hm(xIndex(i),yIndex(i));
        x(i) = xm(xIndex(i),yIndex(i));
        y(i) = ym(xIndex(i),yIndex(i));
    end
    
    pointCloud = [x' y' z'];
    
    save('localTerrainPointCloudRed.mat','pointCloud');
    
end