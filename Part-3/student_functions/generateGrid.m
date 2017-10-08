%% function generate grid
%
% Author Kuan Chun Hwang

function map = generateGrid(mapDim,cellDim,startPos,goalPos)

% Find number of cells required
cell = pos2cellAbs(mapDim,cellDim);

% Input cell colours to white
map = ones(round(cell));

% Find starting position's cell position
startPos_cell = pos2cell(startPos,cellDim,mapDim);
startPos_cellY =round(startPos_cell(1));
startPos_cellX =round(startPos_cell(2));

% Find goal positions's cell position
goalPos_cell = pos2cell(goalPos,cellDim,mapDim);
goalPos_cellY = round(goalPos_cell(1));
goalPos_cellX = round(goalPos_cell(2));

% Make starting colour blue and goal colour red
map(startPos_cellX,startPos_cellY) = 2;
map(goalPos_cellX,goalPos_cellY) = 3;

end