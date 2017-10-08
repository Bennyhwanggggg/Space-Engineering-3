function showMap(fNum,fName,X,Y,map,cMap)
 
    f = figure(fNum);
    colormap(cMap);
    
    image(X,Y,map);
    
    grid on; grid minor; axis image;
    
    xlabel('x [m]')
    ylabel('y [m]')
    
    view(-90,-90);
    
    f.Name = fName;
    f.Units = 'centimeters';
    %f.Position = [0 0 20 20];
    
    drawnow;
    
end