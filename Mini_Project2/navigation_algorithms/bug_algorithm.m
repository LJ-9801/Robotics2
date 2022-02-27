close all; clear all;
% define room
roomspec;
% show room
fignum=1;
h_fig=roomshow(colobj,fignum);
axis('square');
% define robot
rL=.4;rW=.2;rz=.1;
robot=collisionCylinder(rW, rz);
% show robot
q0=[1;2;0];
q(:,1)=q0;


qgoal = [9;5;0];

% ======================================= %
% Bug 0 Algorithm                         %
% Go toward qgoal                         %
% When distance to wpk <= some distance   %
%   record q                              %
%   go around obstacle                    %
%   if q is reach again                   %
%       find min distance around obstacle %
%       go to that position               %
%       go to goal                        %
% ======================================= %


robotshow(robot,q0);

mindist = 0.2;
flag = 0;

k = 1;
while norm(q(1:2,k) - qgoal(1:2)) > 0.5
    [isInt,dist,wp]=colcheck(robot,q(:,k),colobj);
    [rho, ind] = min(dist);
    wparr = [wp{ind}(1,2); wp{ind}(2,2)];

     %vector that connects current state to the witnesspoint nearest to the
     %robot
    vec = [wparr - q(1:2,k); 0];
    vecgoal = qgoal - q(:,k);
    if norm(q(:,k)-qgoal) >= 0.5 && rho >= mindist
        dq = vecgoal./40;
        dq(3) = 0;
        q(:,k+1) = q(:, k)+dq;
        figure(1)
        robotshow(robot,q(:,k+1));
        k=k+1;
    elseif norm(q(:,k)-qgoal) >= 0.5 && rho <= mindist
        [isInt,dist,wp]=colcheck(robot,q(:,k),colobj);
        [rho, ind] = min(dist);
        wparr = [wp{ind}(1,2); wp{ind}(2,2)];
        dq = (vec'*rotz(90))';
        q(:,k+1) = q(:, k)+dq;
        figure(1)
        robotshow(robot,q(:,k+1));
        k=k+1;
    end
    
    
end






