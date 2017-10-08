function cell = pos2cellAbs(mapDim,cellDim)

    xm = mapDim(1);
    ym = mapDim(2);
    cellWidth = cellDim(1);
    cellHeight = cellDim(2);
    
    xg = (xm+cellWidth)/cellWidth;
    yg = (ym+cellHeight)/cellHeight;
    
    cell = [xg,yg];

end