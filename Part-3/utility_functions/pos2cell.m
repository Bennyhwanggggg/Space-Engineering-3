function cell = pos2cell(posCoords,cellDim,mapDim)

    xm = posCoords(1);
    ym = posCoords(2);
    cellWidth = cellDim(1);
    cellHeight = cellDim(2);
    mapWidth = mapDim(1);
    mapHeight = mapDim(2);
    
    xm = xm + mapWidth/2;
    ym = ym + mapHeight/2;
    
    xg = (xm+cellWidth)/cellWidth;
    yg = (ym+cellHeight)/cellHeight;
    
    cell = [xg,yg];

end