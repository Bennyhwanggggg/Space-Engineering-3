%% function roughness hazard evaluation
%
% Author Kuan Chun Hwang

function scoreRoughness = roughnessHazardEval(points,n,p,ROVER_CLEARANCE_HEIGHT)

for i = 1:length(points(:,1))
    r(i) = (points(i,:)-p)*n;
end

std_set = std(r);

scoreRoughness = 255*min(1,(3*std_set)/ROVER_CLEARANCE_HEIGHT);

end