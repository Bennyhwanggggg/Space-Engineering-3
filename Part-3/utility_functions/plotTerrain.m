function plotTerrain(fNum,fName,startPos,goalPos)

load('localTerrain.mat');
hold on
%% Terrain Points
% Plot the triangles represented by all the points. We'll add lighting too
% to make the detail easier to see.
% figure(1);
% trisurf(delaunay(x, y), x, y, h);
% colormap jet;                    % Default color map.
% set(gca, 'Position', [0 0 1 1]); % Fill the figure window.
% axis equal vis3d;                % Set aspect ratio.
% shading interp;                  % Interpolate color across faces.
% camlight left;                   % Add a light over to the left somewhere.
% lighting gouraud;                % Use decent lighting.

%% Meshed Terrain
% Plot as a mesh. Note that this is zoomed in on the median off all of the
% generated points in order to capture the detailed middle instead of the
% relatively uneventful edges.
f = figure(fNum);
f.Name = fName;

surf(xm, ym, hm);
set(gca, 'Position', [0 0 1 1]); % Fill the figure window.
axis equal vis3d;                % Set aspect ratio.
shading interp;                  % Interpolate color across faces.
camlight left;                   % Add a light over to the left somewhere.
lighting gouraud;                % Use decent lighting.
view(168,24);

%Plot start position.
plot3(startPos(1),startPos(2),1.5,...
      'd','MarkerEdgeColor','b','MarkerSize',15,'LineWidth',2);

%Plot goal position.
plot3(goalPos(1),goalPos(2),1.5,...
      'p','MarkerEdgeColor','r','MarkerSize',20,'LineWidth',2);

end