%% Planning for the Orville Mars Rover
%  Author: William Reid(Original) Modified by: Kuan Chun Hwang
%
%% Acknowledgements:
% Some of the planning algorithm material provided in this m-file is
% inspired by the Computational Motion Planning Course offered by CJ Taylor
% from the University of Pennsylvania as part of the Coursera Robotics
% Specialization series of Massive Open Online Courses:
% <https://www.coursera.org/learn/robotics-motion-planning>.
%
% The GESTALT algorithm implementation is based on (Goldberg et al. 2005) Stereo Vision
% and Rover Navigation Software for Planetary Exploration.
%
%% Summary
%
% This m-file provides an introduction to grid-based robot motion planning.
% It is your task to implement a series of traversability map and motion
% planning algorithms for the hypothetical Mars Rover, Orville.
% The topics covered include occupancy grids, configuration spaces, the
% Dijkstra planning algorithm, the A* planning algorithm and a brief look at
% the goodness map construction step of the GESTALT algorithm.
%
%
% This m-file is arranged using cells. You may execute each cell
% individually by clicking on a cell and then pressing the "Run Section"
% button in the Editor banner. Beware that some cells are dependent on variables
% introduced in previous cells.
%
% Let's start things off by clearing the workspace, clearing the command
% window and closing all existing figure windows.

clear; clc; close all;

% Add your functions, utility functions, p_functions and the data set directories
% to the path.
addpath(genpath('p_functions'))
addpath(genpath('student_functions'))
addpath(genpath('utility_functions'))
addpath(genpath('data_sets'))

%% Martian Terrain Map -- There are ZERO functions for you to write in this section.
% Our first task is to define an occupancy grid over which we can plan a
% path for the rover. Each grid cell is to have a size of 20 cm by 20 cm,
% roughly the footprint area of a single Mars Exploration Rover (MER)
% wheel.

%Definition of the start and goal positions.
startPos = [-3.8 4.6];                 %[m]
goalPos = [0.9375 -4.83];              %[m]

%The following plot shows the 10 m by 10 m section of terrain that the
%Orville rover will have to traverse. The start and goal locations of the
%rover are shown with a blue diamond and red star respectively.
plotTerrain(1,'Terrain Map',startPos,goalPos);
xlabel('x [m]');
ylabel('y [m]');
view(56,30)

%The following plot shows a point cloud of the terrain generated from
%previuos exploration of the area by both the Wilbur UAV and the Mars
%Reconnaisance Orbiter (MRO).
figure(2)
plotPointCloud(2,'Terrain Point Cloud');
xlabel('x [m]');
ylabel('y [m]');
view(56,30);

%% Gridded Map -- There is ONE function for you to write in this section.
% Now we will generate a grid that covers the entire terrain. Each grid
% cell will is equal to the footprint area of a single rover wheel:
% approximately 20 cm by 20 cm. It is your task to write a function that
% returns the variable map. map is a matrix with every element of the
% matrix representing a single cell in the occupancy grid.
%
% You shall implement a function generateGrid that returns a
% map matrix that has the correct size to ensure complete coverage of the 20 m by 20 m
% terrain and a cell size of 20 cm by 20 cm. The returned map shall have
% all matrix elements with a value of 1 except for the start and goal
% cells, which will have values of 2 and 3 respectively. These values
% correspond to a enumeration that denotes the type of the cell. This
% enumeration will be used throughout this assignment if not mentioned
% otherwise.
%
% Map grid cell enumeration
% 1 - white   = clear cell
% 2 - blue    = start
% 3 - red     = goal
% 4 - black   = obstacle
% 5 - cyan    = on list
% 6 - orange  = visited
% 7 - yellow  = path

cMapOcc = [1 1 1; 0 0 1; 1 0 0; 0 0 0; 0 1 1; 1 0.5 0; 1 1 0];
cMapTrav = colormap(hot(256));

%Define the dimensions of the local Martian terrain.
MAP_WIDTH  = 10;                %[m] - The distance along the x-axis of the map.
MAP_HEIGHT = 10;                %[m] - The distance along the y-axis of the map.

CELL_WIDTH = 0.2;               %[m] - The width of an individual cell.
CELL_HEIGHT = 0.2;              %[m] - The heigh of an individual cell.

ROVER_FOOTPRINT_RADIUS = 0.5;   %[m] - The radius of the rover's footprint.

mapDim = [MAP_WIDTH MAP_HEIGHT];
cellDim = [CELL_WIDTH CELL_HEIGHT];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map = generateGrid(mapDim,cellDim,startPos,goalPos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mapHeight,mapWidth] = size(map);
numCells = mapWidth*mapHeight;

X = [-MAP_WIDTH/2 MAP_WIDTH/2];
Y = [-MAP_HEIGHT/2 MAP_HEIGHT/2];

showMap(3,'Gridded Map',X,Y,map,cMapOcc);

%% Configuration Space -- There is ONE function for you to write in this section.
% The Orville Rover ground control hazard surveillance team has identified some obvious
% hazard areas from onboard camera data and orbital imagery. They have given the rover
% mobility planning team (you) a list of coordinates in the local terrain
% reference frame that locate the centre of the hazard areas (obstacles).
%
% In this section it is your task to take this list of hazard areas and
% represent them within the map that you generate in the previous section.
% A hazard area or obstacle is indicated in the map using enumeration code
% 4. You will have to write a function generateObstacles that assigns the
% appropriate map cells as 4.
%
% Remember that you are generating the configuration space of the robot and
% not simply the occupancy grid. You must therefore remember to inflate the
% obstacles by a radius equivalent to the radius of the rover's footprint
% (0.5 m).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('hazards.mat')
map = generateObstacles(map,obstacles,cellDim,mapDim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showMap(4,'Configuration Space',X,Y,map,cMapOcc);

%% GESTALT -- There are FIVE functions for you to write in this section.
% Now that we have a configuration space we want to generate a goodness map
% that evaluates the 'free' terrain and returns a traversability map. Each
% cell will be given a traversability score between 0 and 255. This score
% is an estimate of how traversable the associated cell is. The rules for
% estimating this traversability map score are detailed in the Mars Rover
% Planning lecture slides. The lower the cell score, the more traversable
% the cell is.
%
% A point cloud of the local terrain around the Mars Rover is provided in
% the file terrainPointCloud.mat. The GESTALT (Grid-based Estimation of Surface
% Traversability Applied to Local Terrain) algorithm iterates over every
% cell in the map and finds the points within the point cloud that are
% within a radius equivalent to the footprint of the rover. This operation
% is performed by the function getRoverFootprintPoints. You are required to
% implement this function. It is suggested that you look into using the
% rangesearch inbuilt MATLAB function to help implement this function. You
% may find out more about this function by typing 'help rangesearch'.
%
% The next step in GESTALT is to fit a plane to the patch of points
% returned by getRoverFootprintPoints. This operation will be performed by
% the function fitPlane. fitPlane returns a vector, n, that is normal to
% the fitted plane.
%
% With this plane and the patch points, terrain hazard evaluation may be
% conducted. The three hazard evaluation functions to be implemented
% evaluate the incline and roughness of the terrain along with any step
% obstacles that may be present. To evaluate terrain incline the function pitchHazardEval
% finds the angle between the fitted plane's normal vector and the xy plane.
% To evaluate roughness, the function roughnessHazardEval finds the
% standard deviation of the residual of the plane fit. To evaluate step obstacles,
% the function stepHazardEval finds the maximum height difference between points in a cell
% patch. Each of these evaluation functions returns a cell score. The total
% traversability score of cell is taken as the minimum of the incline, roughness and step
% scores.

ROVER_CLEARANCE_HEIGHT = 0.4;   %[m]
ROVER_MAX_PITCH = 25;           %[deg]
TRAVERSABILITY_THRESHOLD = 150; %[-]

mapStep = zeros(mapWidth,mapHeight);
mapRoughness = zeros(mapWidth,mapHeight);
mapPitch = zeros(mapWidth,mapHeight);
mapTraversability = zeros(mapWidth,mapHeight);

for i = 1:numCells
    
    %Get the cell coordinates.
    [yC,xC] = ind2sub(size(map),i);
    
    %Get the position of the cell.
    pos = cell2pos([xC,yC],cellDim,mapDim);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get the indices of the points in the point cloud that are inside the
    %footprint of the rover.
    points = getRoverFootprintPoints(pos,ROVER_FOOTPRINT_RADIUS);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fit a plane to the points.
    [n,p] = fitPlane(points);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Perform a step hazard evaluation over the patch of points.
    mapStep(i) = stepHazardEval(points,ROVER_CLEARANCE_HEIGHT);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Perform a pitch hazard evaluation over the patch of points.
    mapPitch(i) = pitchHazardEval(n,ROVER_MAX_PITCH);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Perform a roughness hazard evaluation over the patch of points.
    mapRoughness(i) = roughnessHazardEval(points,n,p,ROVER_CLEARANCE_HEIGHT);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Use each of the hazard scores to find a traversability score for the
    %cell.
    mapTraversability(i) = max([mapStep(i) mapPitch(i) mapRoughness(i)]);
    
end

%Put the obstacles in the traversability map.
for i = 1:numCells
    if map(i) == 4 || mapTraversability(i) > TRAVERSABILITY_THRESHOLD
        mapTraversability(i) = Inf;
    end
end

%Open a new figure so that the following colormap command doesn't change
%the colormap of the of the figure from the previous section.
showMap(5,'Traversability Map',X,Y,mapTraversability,cMapTrav);


%% Dijkstra's Algorithm -- There are TWO functions for you to write in this section.
% Now that we have a traversability map we can plan a path for the Orville
% Rover. We will first implement Dijkstra's algorithm to perform a breadth
% first search over the traversability map and return a minimum cost path
% between the start and goal nodes.

% Give user option to chose between astar and dijkstra
prompt = {'Dijkstra = 1, A* = 2'};
dlg_title = 'Chose between Dijkstra Algorithm or A*';
answer = inputdlg(prompt,dlg_title);
tr = strcmp('1',answer);
if tr == 1
    dijkstra = true;
else
    dijkstra = false;
end

if dijkstra
    
    %Transform the coordinates of the start and goal nodes into cell indices.
    startCell = pos2cell(fliplr(startPos),cellDim,mapDim);
    goalCell = pos2cell(fliplr(goalPos),cellDim,mapDim);
    startNode = sub2ind(size(map),int64(startCell(1)),int64(startCell(2)));
    goalNode = sub2ind(size(map),int64(goalCell(1)),int64(goalCell(2)));
    
    %Make this flag true to enable visualisation of the planner on every
    %iteration of the planning loop.
    % Give user option to visualise the planner
    prompt = {'Visualise the planner (drawmap): Yes or No?'};
    dlg_title = 'Dijkstra Planner visualisation';
    answer = inputdlg(prompt,dlg_title);
    dr = strcmp('Yes',answer);
    if dr == 1
        drawMapEveryTime = true;
    else
        drawMapEveryTime = false;
    end
    
    % Initialize the cost array
    distanceFromStart = Inf*ones(mapWidth,mapHeight);
    distanceFromStart(startNode) = 0;
    minDist = 0;
    
    %Initialize the exploration at the start node.
    currentNode = startNode;
    
    %Set the parent of the initial node to null.
    parent = zeros(mapWidth,mapHeight);
    
    % Keep track of number of nodes expanded
    numExpanded = 0;
    
    %Initialize the open node array, H.
    H = Inf*ones(mapWidth,mapHeight); %try no inf
    H(startNode) = 1;
    
    %Initialize the closed node array, W.
    W = zeros(mapWidth,mapHeight);
    
    fName = 'Dijkstra ''s Algorithm';
    
    % Main planning loop
    while (min(W(:)) < 1 && currentNode ~= goalNode && ~isinf(minDist))
        
        % Draw current map. Make drawMapEveryTime = true if you want to see how the
        % nodes are expanded on the grid.
        if (drawMapEveryTime)
            showMap(6,fName,X,Y,map,cMapOcc)
        end
        
        % Find the node with the minimum distance
        tempMat = H.*distanceFromStart;
        [minDist,currentNode] = min(tempMat(:));
        
        % Update map
        map(currentNode) = 6; % mark current node as visited
        
        % Compute row, column coordinates of current node
        [i, j] = ind2sub(size(distanceFromStart), currentNode);
        
        %Iterate numExpanded
        numExpanded = numExpanded+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Visit each neighbor of the current node and update the map, distances
        % and parent tables appropriately.
        conn = 8; %The neighbour-connectedness, can be either 4 or 8 connected.
        
        %getNeighbours finds the valid neighbours surrounding the cell at
        %index (i,j). map is the occupancy grid, mapHeight and mapWidth are
        %the dimensions of the map and conn is the neighbour connectedness
        %(either 4 or 8). The outputs of this function are neighbours,
        %which is a vector the cell linear indices of the neighbouring
        %nodes. neighbourInd is a vector of the neighbour indices relative
        %to the centre node. This vector can be used to help identify which
        %neighbouring nodes exist within the Dijkstra and A* algorithms.
        %You can output neighbourInd as [] if you do not want to use it.
        [neighbours,neighbourInd] = getNeighbours(i,j,map,mapHeight,mapWidth,conn);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [distanceFromStart,parent,H,map] = ...
            dijkstraBody(currentNode,distanceFromStart,parent,H,map,mapTraversability,neighbours,neighbourInd);
        
        
        minDist = distanceFromStart(currentNode);
        W(currentNode) = 1;
        H(currentNode) = Inf;
        
    end
    
    if (isinf(distanceFromStart(goalNode)))
        route = [];
        cost = Inf;
    else
        route = goalNode;
        cost = 0;
        
        while (parent(route(1)) ~= 0)
            cost = cost + distanceFromStart(parent(route(1)));
            route = [parent(route(1)), route];
        end
        
        % Snippet of code used to visualize the map and the path
        for k = 1:length(route)
            map(route(k)) = 7;
            showMap(6,fName,X,Y,map,cMapOcc);
        end
        
        fprintf('Dijkstra: Cost = %d\n',cost);
    end
    
    map(startNode) = 2;
    map(goalNode) = 3;
    
    showMap(6,fName,X,Y,map,cMapOcc);
    
    
    %% A* Algorithm -- There is ONE function to write in this section.
    
else
    
    %Transform the coordinates of the start and goal nodes into cell indices.
    startCell = pos2cell(fliplr(startPos),cellDim,mapDim);
    goalCell = pos2cell(fliplr(goalPos),cellDim,mapDim);
    startNode = sub2ind(size(map),int64(startCell(1)),int64(startCell(2)));
    goalNode = sub2ind(size(map),int64(goalCell(1)),int64(goalCell(2)));
    
    %Make this flag true to enable visualisation of the planner on every
    %iteration of the planning loop.
    % Give user option to visualise the planner
    prompt = {'Visualise the planner (drawmap): Yes or No?'};
    dlg_title = 'A* Planner visualisation';
    answer = inputdlg(prompt,dlg_title);
    dr = strcmp('Yes',answer);
    if dr == 1
        drawMapEveryTime = true;
    else
        drawMapEveryTime = false;
    end
    
    % Initialize the cost array
    minDist = 0;
    
    %Initialize the exploration at the start node.
    currentNode = startNode;
    
    %Set the parent of the initial node to null.
    parent = zeros(mapWidth,mapHeight);
    
    % Keep track of number of nodes expanded
    numExpanded = 0;
    
    %Initialize the open node array, H.
    H = Inf*ones(mapWidth,mapHeight);
    H(startNode) = 1;
    
    %Initialize the closed node array, W.
    W = zeros(mapWidth,mapHeight);
    
    % Evaluate Heuristic function, H, for each grid cell.
    [XM, YM] = meshgrid (1:mapWidth, 1:mapHeight);
    xd = goalCell(2);
    yd = goalCell(1);
    
    %Manhattan distance heuristic
    %heuristic = abs(XM - xd) + abs(YM - yd);
    %Euclidean distance heuristic
    heuristic = sqrt((XM-xd).^2 + (YM-yd).^2);
    
    fName = 'A* Algorithm';
    
    % Initialize cost arrays
    f = Inf(mapWidth,mapHeight);
    g = Inf(mapWidth,mapHeight);
    g(startNode) = 0;
    f(startNode) = heuristic(startNode);
    
    % Main Loop
    while (min(W(:)) < 1 && currentNode ~= goalNode && ~isinf(minDist))
        
        % Draw current map. Make drawMapEveryTime = true if you want to see how the
        % nodes are expanded on the grid.
        if (drawMapEveryTime)
            showMap(6,fName,X,Y,map,cMapOcc);
        end
        
        % Find the node with the minimum distance
        tempMat = H.*f;
        [~,currentNode] = min(tempMat(:));
        
        % Update map
        map(currentNode) = 6; % mark current node as visited
        
        % Compute row, column coordinates of current node
        [i, j] = ind2sub(size(g), currentNode);
        
        %Iterate numExpanded
        numExpanded = numExpanded+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        conn = 8; %The neighbour-connectedness, can be either 4 or 8 connected.
        [neighbours,neighbourInd] = getNeighbours(i,j,map,mapHeight,mapWidth,conn);
        

        [g,f,parent,H,map] = ...
            aStarBody(currentNode,g,f,parent,H,heuristic,map,mapTraversability,neighbours,neighbourInd);
        

        W(currentNode) = 1;
        H(currentNode) = Inf;
    end
    
    if (isinf(heuristic(goalNode)))
        route = [];
        cost = Inf;
    else
        route = goalNode;
        cost = 0;
        
        while (parent(route(1)) ~= 0)
            cost = cost + g(parent(route(1)));
            route = [parent(route(1)), route];
        end
        
        %Visualize the path.
        for k = 1:length(route)
            map(route(k)) = 7;
            showMap(6,fName,X,Y,map,cMapOcc);
        end
        
        fprintf('Cost = %d\n',cost);
    end
    
    map(startNode) = 2;
    map(goalNode) = 3;
    
    showMap(6,fName,X,Y,map,cMapOcc);
    
end

%% Plot Path on Terrain
figure(1)

load('localTerrainPointCloud');

for i = 1:length(route)
    [yC,xC] = ind2sub(size(map),route(i));
    pos = cell2pos([xC yC],cellDim,mapDim);
    z = pointCloud(knnsearch(pointCloud,[pos 0],'K',1),3)+0.2;
    plot3(pos(1),pos(2),z,'r*-');
end

%% Extension Geological Survey Expedition
% Give user option to do geological expedition
prompt = {'Extension Geological Survey Expedition: Yes or No?'};
dlg_title = 'Geological Survey Expedition';
answer = inputdlg(prompt,dlg_title);
ext = strcmp('Yes',answer);
if ext == 1
    extension = true;
else
    extension = false;
end

% For geological survey expedition the A* algorithm has been implemented
if extension
    % Load the site of interest file
    load('sitesOfInterest.mat')
    
    % Retrieve the map with obstacles, traversability map and basic path
    % planner for the Orville rover
    map2 = map;
    
    % Clear out the old map except for the obstacles
    map2(map2==2) = 1;
    map2(map2==3) = 1;
    map2(map2==5) = 1;
    map2(map2==6) = 1;
    map2(map2==7) = 1;
    
    % Store sites of interest and display it on the map with color blue
    sitesofInterest = [x y'];
    [map2, siteOfinterestPos] = plotSiteOfInterest(sitesofInterest,map2,cellDim,mapDim);
    
    %Make this flag true to enable visualisation of the planner on every
    %iteration of the planning loop.
    % Give user option to visualise the planner
    prompt = {'Visualise the planner (drawmap): Yes or No?'};
    dlg_title = 'Extension Planner visualisation';
    answer = inputdlg(prompt,dlg_title);
    dr = strcmp('Yes',answer);
    if dr == 1
        drawMapEveryTime = true;
    else
        drawMapEveryTime = false;
    end
    
    % Initialize route set
    route_set = 1;
    
    % The double for loop goes through every possible combination as b
    % always have to be larger than a which means no repeated paths are
    % calculated and the start position and goal position would never be
    % the same point. Only one direction of the path is required since if
    % the start position and goal position are reversed, the path would
    % still be the same
    for a = 1:length(siteOfinterestPos(:,1))
        for b = 1:length(siteOfinterestPos(:,1))
            if a<b
                
                % Clear out the old coloured cells
                map2(map2==5) = 1;
                map2(map2==6) = 1;
                
                %Transform the coordinates of the start and goal nodes into cell indices.
                startCell = pos2cell(fliplr(siteOfinterestPos(a,:)),cellDim,mapDim);
                goalCell = pos2cell(fliplr(siteOfinterestPos(b,:)),cellDim,mapDim);
                startNode = sub2ind(size(map2),int64(startCell(1)),int64(startCell(2)));
                goalNode = sub2ind(size(map2),int64(goalCell(1)),int64(goalCell(2)));
                
                % Initialize the cost array
                minDist = 0;
                
                %Initialize the exploration at the start node.
                currentNode = startNode;
                
                %Set the parent of the initial node to null.
                parent = zeros(mapWidth,mapHeight);
                
                % Keep track of number of nodes expanded
                numExpanded = 0;
                
                %Initialize the open node array, H.
                H = Inf*ones(mapWidth,mapHeight);
                H(startNode) = 1;
                
                %Initialize the closed node array, W.
                W = zeros(mapWidth,mapHeight);
                
                % Evaluate Heuristic function, H, for each grid cell.
                [XM, YM] = meshgrid (1:mapWidth, 1:mapHeight);
                xd = goalCell(2);
                yd = goalCell(1);
                
                %Manhattan distance heuristic
                %heuristic = abs(XM - xd) + abs(YM - yd);
                %Euclidean distance heuristic
                heuristic = sqrt((XM-xd).^2 + (YM-yd).^2);
                
                fName = 'A* Algorithm Sites of Interest';
                
                % Initialize cost arrays
                f = Inf(mapWidth,mapHeight);
                g = Inf(mapWidth,mapHeight);
                g(startNode) = 0;
                f(startNode) = heuristic(startNode);
                
                % Main Loop
                while (min(W(:)) < 1 && currentNode ~= goalNode && ~isinf(minDist))
                    
                    % Draw current map. Make drawMapEveryTime = true if you want to see how the
                    % nodes are expanded on the grid.
                    if (drawMapEveryTime)
                        showMap(7,fName,X,Y,map2,cMapOcc);
                    end
                    
                    % Find the node with the minimum distance
                    tempMat = H.*f;
                    [~,currentNode] = min(tempMat(:));
                    
                    % Update map
                    map2(currentNode) = 6; % mark current node as visited
                    
                    % Compute row, column coordinates of current node
                    [i, j] = ind2sub(size(g), currentNode);
                    
                    %Iterate numExpanded
                    numExpanded = numExpanded+1;
                    
                    conn = 8; %The neighbour-connectedness, can be either 4 or 8 connected.
                    [neighbours,neighbourInd] = getNeighbours(i,j,map2,mapHeight,mapWidth,conn);
                    
                    [g,f,parent,H,map2] = ...
                        aStarBody(currentNode,g,f,parent,H,heuristic,map2,mapTraversability,neighbours,neighbourInd);
                    
                    W(currentNode) = 1;
                    H(currentNode) = Inf;
                end
                
                if (isinf(heuristic(goalNode)))
                    route = [];
                    cost = Inf;
                else
                    route = goalNode;
                    cost = 0;
                    
                    while (parent(route(1)) ~= 0)
                        cost = cost + g(parent(route(1)));
                        route = [parent(route(1)), route];
                    end
                    route_num{route_set} = route;
                    
                    route_set = route_set + 1;
                    
                    % Print the cost of the path
                    fprintf('Cost to site of interest = %d\n',cost);
                    
                    figure(8)
                    plotTerrain(8,'Terrain Map',siteOfinterestPos(a,:),siteOfinterestPos(b,:));
                    xlabel('x [m]');
                    ylabel('y [m]');
                    view(56,30)
                    
                    
                end
            end
            
        end
        
    end
    
    % Show the path on 2D figure if drawMapEveryTime is true
    if (drawMapEveryTime)
        
        %Visualize the path.
        for rl = 1: length(route_num(1,:))
            
            % Extract the route
            path = route_num{rl};
            
            for k = 1:length(path)
                
                % Show yellow colour on path
                map2(path(k)) = 7;
                
                % Show the sites of interest on the map
                [map2, siteOfinterestPos] = plotSiteOfInterest(sitesofInterest,map2,cellDim,mapDim);
                showMap(7,fName,X,Y,map2,cMapOcc);
                
            end
        end
    end
    
    % Load terrain point cloud file
    load('localTerrainPointCloud');
    
    % Plot the set of paths on to the terrain
    for rl = 1: length(route_num(1,:))
        path = route_num{rl};
        for k = 1:length(path)
            figure(8)
            [yC,xC] = ind2sub(size(map2),path(k));
            pos = cell2pos([xC yC],cellDim,mapDim);
            z = pointCloud(knnsearch(pointCloud,[pos 0],'K',1),3)+0.2;
            plot3(pos(1),pos(2),z,'r*-');
        end
    end
end

