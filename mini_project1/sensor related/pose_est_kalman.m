%
% pose_est_kalman.m
%
% Given y and robot kinematics, use Extended Kalman Filter (EKF) to
% estimate the robot position and orientation 
%
%function qhat=pose_est_kalman(y,pL,pB,N_scan,wcov,vcov)

function [qhatnext,Pnext,Snext]=... 
    pose_est_kalman(qhat,u,y,ts,pL,pB,N_scan,wcov,vcov,P,Sigma,robot,colobj)

% use range sensors first
N_range=size(pL,2);N_bearing=size(pB,2);N_odo=2;

% propagate the state via dead reckoning to the next state using current
% estimate 

[qbarnext,utrue]=wmr(qhat,u,ts,(zeros(size(wcov))));
qbarnext=qbarnext';

yind=[(1:N_range) (N_range+N_odo+1:N_range+N_odo+N_bearing)];
A=[1 0 -u(1)*sin(qhat(3));0 1 u(1)*cos(qhat(3)); 0 0 1];
B=[cos(qhat(3)) 0 ; sin(qhat(3)) 0 ; 0 1 ];
V=diag(vcov(yind));
W=diag(wcov);

C=[gradient_range(pL,qbarnext);gradient_bearing(pB,qbarnext)];

K=P*C'*inv(V);
ybar=h(qbarnext,pL,pB);

qhatnext=qbarnext+K*(y(yind)-ybar);
Pnext=inv(inv(Sigma)+C'*inv(V)*C);
Snext=A*Pnext*A'+B*W*B';

end

% range and bearing sensor output map

function y=h(q,pL,pB)
    pR=[q(1:2);0];
    y_range=vecnorm(pR-pL)';
    y_bearing=(atan2(pB(2,:)-q(2),pB(1,:)-q(1))-q(3))';
    y=[y_range;y_bearing];
end

% gradient of range sensor model

function grad_h=gradient_range(pL,q)
    pR=[q(1:2);0];
    grad_h=(-(pL-pR)./vecnorm(pL-pR))';    
end

% gradient of bearing sensor model

function grad_h=gradient_bearing(pB,q)
    pR=[q(1:2);0];
    grad_h=[([pB(2,:)'-pR(2) -(pB(1,:)'-pR(1))]'./vecnorm(pB(1:2,:)-pR(1:2)).^2)' ...
        -ones(size(pB,2),1)];
end