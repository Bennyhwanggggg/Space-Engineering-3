%% function find best fit plane
%
% Author Kuan Chun Hwang

function [n,p] = fitPlane(points)

% Find average points
p = [mean(points(:,1)),mean(points(:,2)),mean(points(:,3))];

% Calculate R matrix
for i = 1:length(points(:,1))
    R(i,:) = points(i,:)-p;
    
end

% Find eigenvalues and eigenvectors of the matrix R'*R
[V,D] = eig(R'*R);

% Fuind normal vector corresponding to the best fit plane and the unit
% normal vector
normal = V(:,1);
n = normal/norm(normal);

end