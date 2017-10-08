%% function step hazard evaluation
%
% Author Kuan Chun Hwang

function scoreStep = stepHazardEval(points,ROVER_CLEARANCE_HEIGHT)

% Find maximum height difference
maxHeightDifference = peak2peak(points(:,3));

% see if the step is above clearance height or else assign the hazard
% evaluation score
if maxHeightDifference <ROVER_CLEARANCE_HEIGHT/3
    scoreStep = 0;
else
    scoreStep = 255*min(1,maxHeightDifference/ROVER_CLEARANCE_HEIGHT);
end

end