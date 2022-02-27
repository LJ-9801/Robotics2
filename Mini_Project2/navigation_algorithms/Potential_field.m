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
q0=[.2;2;0];
q(:,1)=q0;



qgoal = [9;5;0];

Kp = 0.05;
eta = 0.05;
alpha = 0.5;
rho0 = 0.5;
delta = 0.01;


robotshow(robot,q0);



%===================================================
%q is the state of the robot
%define U(qk) = Uatt + Urep
%Uatt = 0.5*Kp||q-qgoal||^2
%rho(q) is the distance to the minimum object distance
%
% while the robot has not reached the goal
k = 1;
while norm(q(1:2,k) - qgoal(1:2)) > 0.5

    % attractive force
    Uatt = gradU_att(q(:,k), qgoal, Kp);
    % repulsive force
    gradUrep = gradU_rep(q(:,k), robot, colobj, rho0, eta);
    %gradient descent with infinite wall
    dq = -Uatt-alpha*gradUrep';
    %q(:,k+1) = q(:,k)+[dq;0];
    q(:,k+1) = q(:,k)+[dq;0]+ (rand-0.5)*2*[0.09;0.09;0];
    [isInt,dist,wp]=colcheck(robot,q(:,k),colobj);
    [rho, ind] = min(dist);
    wparr(:,k) = [wp{ind}(1,2); wp{ind}(2,2)];
    plot(wparr(1,k), wparr(2,k),'kx', 'LineWidth', 10);
    figure(1)
    robotshow(robot,q(:,k+1));
    k=k+1;
end



function Uatt = gradU_att(qk, qgoal, Kp)
    Uatt = Kp*(qk(1:2)-qgoal(1:2));
end



function gradUrep = gradU_rep(q, robot, colobj, rho0, eta)
    e = 0.01;
    [isInt,dist,wp]=colcheck(robot,q,colobj);
    [rho, ind] = min(dist);

    [isInt,dist1,wp]=colcheck(robot,q+e*[1;0;0],colobj);
    [rho1, ind] = min(dist1);

    [isInt,dist2,wp]=colcheck(robot,q+e*[0;1;0],colobj);
    [rho2, ind] = min(dist2);

    gradrho = [(rho1-rho)/e (rho2-rho)/e];

    gradUrep = [0 0];
    if rho<rho0
        gradUrep = (-eta/rho^2)*(1/rho-1/rho0)*gradrho;
    end

end