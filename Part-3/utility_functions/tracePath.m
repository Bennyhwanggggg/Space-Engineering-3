function path = tracePath(costMap,currNode)

    [i,j] = ind2sub(size(costMap),currNode);
    [h,w] = size(costMap);

    while costMap(currNode) > 0
    
        if i-1 > 1 && i-1 <= w && j > 1 && j <= h
            if costMap(i-1,j) == costMap(currNode) - 1
                currNode = sub2ind(size(costMap),i-1,j);
                path = [path currNode];
            end
        elseif i > 1 && i <= w && j-1 > 1 && j-1 <= h
            if costMap(i-1,j) == costMap(currNode) - 1
                currNode = sub2ind(size(costMap),i,j-1);
                path = [path currNode];
            end
        elseif i+1 > 1 && i+1 <= w && j > 1 && j <= h
            if costMap(i-1,j) == costMap(currNode) - 1
                currNode = sub2ind(size(costMap),i+1,j);
                path = [path currNode];
            end
        elseif i > 1 && i <= w && j+1 > 1 && j+1 <= h
            if costMap(i-1,j) == costMap(currNode) - 1
                currNode = sub2ind(size(costMap),i,j+1);
                path = [path currNode];
            end
        else
            %No path found.
            break;
            path = [];
        end
        
    end
    
end