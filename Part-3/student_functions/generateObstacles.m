%% function generate obstacles
%
% Author Kuan Chun Hwang

function map = generateObstacles(map,obstacles,cellDim,mapDim)

% Radius of rover's footprint
radius = 0.5;

% Find the footprint radius interms of cell. Assume the dimension of cells
% are always squares
r_cell = radius./cellDim;

% Create x and y coordinates in terms of map size using meshgrid
[xx, yy] = meshgrid(1:length(map(:,1)),1:length(map(:,2)));

% Use a loop to plot position of each obstacle
for i = 1:size(obstacles,1)
    
    % Find obstacle center cell positions
    obstacles_cell = pos2cell(obstacles(i,1:2),cellDim,mapDim);
    
    % Use circles formula to calculate the obstacles cells with rover
    % footprint radius and the obstacles radius taken into account
    obstacles_r = sqrt((xx-obstacles_cell(1)).^2 + (yy-obstacles_cell(2)).^2)<=(r_cell(1)+obstacles(i,3)./cellDim(1));
    
    % Find these positions in the map matrix
    obstacle_loc = find(obstacles_r==1);
    
    % Change obstacle cells' colour to black
    map(obstacle_loc) = 4;
end

end