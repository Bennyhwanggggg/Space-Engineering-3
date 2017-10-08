%% function 
%
% Author: Kuan Chun Hwang

function [g,f,parent,H,map] = aStarBody(currentNode,g,f,parent,H,heuristic,map,mapTraversability,neighbours,neighbourInd)

% Combine neighbour cells with their index positions
nei_cell = [neighbours neighbourInd];

% Compute the cost of each neighbour cells and compute the cost to come if it is a shorter path with lower cost.
for i = 1: length(neighbours(:,1))
    
    % Diagnol cells have a longer distance which is a ratio of sqrt(2)
    if nei_cell(i,3) == 1 || nei_cell(i,3) == 2 || nei_cell(i,3) == 3 ||nei_cell(i,3) == 4
        cost = sqrt(2)*(mapTraversability(neighbours(i,1),neighbours(i,2))+mapTraversability(currentNode))/2;
    else
        cost = (mapTraversability(neighbours(i,1),neighbours(i,2))+mapTraversability(currentNode))/2;
    end
    
    if  g(neighbours(i,1),neighbours(i,2)) > g(currentNode) + cost
        g(neighbours(i,1),neighbours(i,2)) = g(currentNode) + cost;
         % Assign the cyan color on map for that neighbouring cell for
        % visualisation
        map(neighbours(i,1),neighbours(i,2)) = 5;
        
        % Set open set
        H(neighbours(i,1),neighbours(i,2)) = 1;
        
        
        % Set neighbouring cell's parent to the current node
        parent(neighbours(i,1),neighbours(i,2)) = currentNode;
        
        % Compute f which is the total of cost to come and cost to go
        f(neighbours(i,1),neighbours(i,2)) = g(neighbours(i,1),neighbours(i,2)) + heuristic(neighbours(i,1),neighbours(i,2));
    end
end
end