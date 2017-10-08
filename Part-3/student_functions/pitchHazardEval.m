%% function pitch hazard evaluation
%
% Author Kuan Chun Hwang

function scorePitch = pitchHazardEval(n,ROVER_MAX_PITCH)

% Normal vector to x-y plane
n_xy = [0 0 1];

% Find dot product of the two normals
c_ang = dot(n',n_xy);

% Find the angle between two planes
pitch = acosd((c_ang/(norm(n_xy)*norm(n))));

% See if the pitch is larger than maximum rover's pitch assign the hazard value to 255, else assign hazard
% evaluation value
if pitch > ROVER_MAX_PITCH
    scorePitch = 255;
else
    scorePitch = 255*min(1,pitch/ROVER_MAX_PITCH);
end

end


