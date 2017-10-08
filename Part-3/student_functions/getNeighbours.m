%% function get beighbouring cell indices
%
% Author Kuan Chun Hwang

function [neighbours,neighbourInd] = getNeighbours(i,j,map,mapHeight,mapWidth,conn)

% Initialize row number
a = 1;

% Find the cell that was already used as center and the cells occupied by
% obstacles
[i_u, j_u] = find(map==6);
[i_ob, j_ob] = find(map==4);

% for 8 connected cell go through each neighbouring cell by using i, i+1,
% i-1, j,j+1 and j-1. By considering the map size as a condition as well,
% so that when it is exceeding the map, the cell is not computed. Each
% neighbour cells are assigned a index position as well
if conn == 8
    for y_ind = 1:mapWidth
        if y_ind == i
            if j-1 >= 1
                if map(i,j-1) ~= 4
                    neighbours(a,:) = [i j-1];
                    neighbourInd(a,:) = 6;
                    a = a + 1;
                end
            end
            if mapHeight >= j+1
                if map(i,j+1) ~= 4
                    neighbours(a,:) = [i j+1];
                    neighbourInd(a,:) = 8;
                    a = a + 1;
                end
            end
        end
        if y_ind == i-1 && i-1 >= 1
            
            if map(i-1,j) ~= 4
                neighbours(a,:) = [i-1 j];
                neighbourInd(a,:) = 5;
                a = a + 1;
            end
            
            if j-1 >= 1
                if map(i-1,j-1) ~= 4
                    neighbours(a,:) = [i-1 j-1];
                    neighbourInd(a,:) = 2;
                    a = a + 1;
                end
            end
            if mapHeight >= j+1
                if map(i-1,j+1) ~= 4
                    neighbours(a,:) = [i-1 j+1];
                    neighbourInd(a,:) = 3;
                    a = a + 1;
                end
            end
        end
        if y_ind == i+1 && mapWidth >= i+1
            if map(i+1,j) ~= 4
                neighbours(a,:) = [i+1 j];
                neighbourInd(a,:) = 7;
                a = a + 1;
            end
            if  j-1 >= 1
                if map(i+1,j-1) ~= 4
                    neighbours(a,:) = [i+1 j-1];
                    neighbourInd(a,:) = 1;
                    a = a + 1;
                end
            end
            if mapHeight >= j+1
                if map(i+1,j+1) ~= 4
                    neighbours(a,:) = [i+1 j+1];
                    neighbourInd(a,:) = 4;
                    a = a + 1;
                end
            end
        end
    end
    
% For 4 connected, only compute the right, left, up and down cells    
elseif conn == 4
    for y_ind = 1:mapWidth
        if y_ind == i
            if j-1 >= 1
                if map(i,j-1) ~= 4
                    neighbours(a,:) = [i j-1];
                    neighbourInd(a,:) = 6;
                    a = a + 1;
                end
            end
            if mapHeight >= j+1
                if map(i,j+1) ~= 4
                    neighbours(a,:) = [i j+1];
                    neighbourInd(a,:) = 8;
                    a = a + 1;
                end
            end
        end
        if y_ind == i-1 && i-1 >= 1
            
            if map(i-1,j) ~= 4
                neighbours(a,:) = [i-1 j];
                neighbourInd(a,:) = 5;
                a = a + 1;
            end
        end
        
        if y_ind == i+1 && mapWidth >= i+1
            if map(i+1,j) ~= 4
                neighbours(a,:) = [i+1 j];
                neighbourInd(a,:) = 7;
                a = a + 1;
            end
        end
        
    end
end

% Combine the neighbour cells with their index position and sort it in
% descending order based on index positions
neigh_combined = [neighbours neighbourInd];
sorted_neighbours = sortrows(neigh_combined,-3);
neighbours = sorted_neighbours(:,1:2);
neighbourInd = sorted_neighbours(:,3);

% Remove the cells that were already visitied or occupied obstacles
repPathCell = find(ismember(neighbours(:,1:2),[i_u j_u],'rows'));
obstacleCell = find(ismember(neighbours(:,1:2),[i_ob j_ob],'rows'));
if isempty(repPathCell) == 0
    neighbours(repPathCell,:) = [];
    neighbourInd(repPathCell,:) = [];
end
if isempty(obstacleCell) == 0
    neighbours(obstacleCell,:) = [];
    neighbourInd(obstacleCell,:) = [];
end


end