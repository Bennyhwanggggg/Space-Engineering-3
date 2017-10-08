%% Plot orbit transfer
% This function plots orbit transfer using universal conic section and also
% the initial and target orbit using ECI positions
%
% Author: Kuan Chun Hwang

function plotOrbitTransfer(pos_eci_initial,pos_eci_final,Optimal_dE1,r0,v0,Optimal_dE2,Optimal_theta1,Optimal_psi1,Optimal_dV1,Optimal_dV2,Optimal_theta2,Optimal_psi2)

%% Plot initial orbit and final orbit using known orbital parameters
figure
plotEarth
hold on
initial_orbit = plot3(pos_eci_initial(:, :, 1), pos_eci_initial(:, :, 2), pos_eci_initial(:, :, 3),'y');
hold on
Target_orbit = scatter3(pos_eci_final(:, :, 1), pos_eci_final(:, :, 2), pos_eci_final(:, :, 3),'.c');
hold on

%% Initial Coast
% Find initial orbit before first burn
for angle = 1:Optimal_dE1
    
    % Use universal conic section method to find all the velocity and
    % position during coast
    [ r1(:,angle),v1(:,angle) ] = UniConic( r0,v0, deg2rad(angle) );
    
end

%% First burn to tranfer orbit
% Add the optimal deltaV1 value to initiate first burn
dV1_eci = deltaV(Optimal_dV1,v1(:,end),r1(:,end),deg2rad(Optimal_theta1),deg2rad(Optimal_psi1));
% Find the new velocity after first burn
v2 = v1(:,end) + dV1_eci;
% The position remains the same assuming instanetaneous burn
r2=r1(:,end);

% Mark this position on plot and show the velcoity vector applied at this
% burn (Vector is scaled by 2 for visibility)
FirstBurnPos = plot3(r2(1),r2(2),r2(3),'*','MarkerSize',15);
hold on
FirstBurnVelVector = quiver3(r2(1),r2(2),r2(3),2*norm(dV1_eci)*dV1_eci(1),5*norm(dV1_eci)*dV1_eci(2),5*norm(dV1_eci)*dV1_eci(3),'r');
hold on

% Find transfer orbit positions and velocity vectors
for angle=1:Optimal_dE2
    
    % % Use universal conic section method to find all the velocity and
    % position during transfer
    [ r3(:,angle),v3(:,angle) ] = UniConic(r2,v2,deg2rad(angle));
    
end

%% Second burn to transit into final orbit
% Add the optimal deltaV2 value to initiate second burn
dV2_eci = deltaV(Optimal_dV2,v3(:,end),r3(:,end),deg2rad(Optimal_theta2),deg2rad(Optimal_psi2));
% Find the new velocity after second burn
v4 = v3(:,end) + dV2_eci;
% The position remains the same assuming instanetaneous burn
r4 = r3(:,end);

% Mark this position on plot and show the velcoity vector applied at this
% burn (Vector is scaled by 2 for visibility)
SecondBurnPos = plot3(r4(1),r4(2),r4(3),'o','MarkerSize',15);
hold on
SecondBurnVelVector = quiver3(r4(1),r4(2),r4(3),2*norm(dV2_eci)*dV2_eci(1),5*norm(dV2_eci)*dV2_eci(2),5*norm(dV2_eci)*dV2_eci(3),'r');
hold on

% Find final orbit after second burn
for angle=1:(360*2)
    
    % % Use universal conic section method to find all the velocity and
    % position in final orbit
    [ r5(:,angle),v5(:,angle) ] = UniConic(r4,v4,deg2rad(angle));
    
end

%% Orbit plotting
% Plot initial coasting orbit
Initial_Orbit = plot3(r1(1,:),r1(2,:),r1(3,:),'b');
hold on
% Plot transfer orbit
Tranfer_orbit=plot3(r3(1,1:end),r3(2,1:end),r3(3,1:end),'k');
hold on
% Plot final orbit achieved
Final_orbit=plot3(r5(1,:),r5(2,:),r5(3,:),'r','LineWidth',3);

% Graph properties
legend([initial_orbit,Initial_Orbit,Tranfer_orbit,Target_orbit,Final_orbit,FirstBurnPos,SecondBurnPos,FirstBurnVelVector,SecondBurnVelVector],'Satellie Parking Orbit Coast','Initial Satellite Parking Orbit','Transfer Orbit','Target Orbit','Final Geostationary Orbit achieved','1st Burn Position','2nd Burn Position','1st Burn Velocity Vector', '2nd Burn Velocity Vector')
title('Orbit Transfer')
view(142.5,30)
end