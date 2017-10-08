function pos = cell2pos(gridCoords,cellDim,mapDim)

    xg = gridCoords(1); yg = gridCoords(2);
    cellWidth = cellDim(1); cellHeight = cellDim(2);
    mapWidth = mapDim(1); mapHeight = mapDim(2);
    
    xm = xg*cellWidth-cellWidth;
    ym = yg*cellHeight-cellHeight;
    
    xm = xm - mapWidth/2;
    ym = ym - mapHeight/2;
    
    pos = [xm ym];

end