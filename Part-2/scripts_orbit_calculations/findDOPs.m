%% function find DOPs
%
% Author Kuan Chun Hwang

function DOPs = findDOPs(H_set)

         V = inv(H_set'*H_set);
         DOPvs = diag(V);
         Vx = DOPvs(1);
         Vy = DOPvs(2);
         Vz = DOPvs(3);
         Vt = DOPvs(4);
         GDOP = sqrt(Vx+Vy+Vz+Vt);
         PDOP = sqrt(Vx+Vy+Vz);
         HDOP = sqrt(Vx+Vy);
         VDOP = sqrt(Vz);
         TDOP = sqrt(Vt);
         
         DOPs = [GDOP;PDOP;HDOP;VDOP;TDOP];
end