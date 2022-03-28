%
% probability road map (PRM) method
% 
% to run: just type in 
% 
% >> prm_move
% 

clear all;close all;

% define room
roomspec;
% generate occupany map
resolution = 100;
room_map = colobj2binary_map(colobj, resolution);

% define robot
rL=.4;rW=.2;rz=.1;
% for formation
rW=.6;
robot=collisionCylinder(rW,rz);

% find configuration space
cspace_map = copy(room_map);
inflate(cspace_map,rW);

% generate PRM
room_prm=mobileRobotPRM(cspace_map,10);

% initial condition
q0=[rW+.1;2;0];

more_qf=1;

rng(1000);
while more_qf>0
    % check if random qf should be used
    qf_flag=1;
    while qf_flag>0
        qf=[5*rand(1,2)+[5 5] 2*pi*(rand-.5)]';
        [isInt,dist,wp]=colcheck(robot,qf,colobj);
        if (max(isnan(dist))==0)&&(min(dist)>rW);qf_flag=0;end
    end
    disp(sprintf('qf = [%0.3g, %0.3g]', qf(1:2)));
    % generate feasible path
    path = findpath(room_prm,q0(1:2)',qf(1:2)');

    figure(1);
    show(room_prm); hold on
    roomshow(colobj, 1);
    axis('square');
    path_L=sum(vecnorm(diff(path)'));
    for i=1:length(path)
        if i==1
            theta(i)=q0(3);
        else
            path_L_i=sum(vecnorm(diff(path(1:i,:))'));
            theta(i)=(qf(3)-q0(3))*path_L_i/path_L+q0(3);
        end
        robotshow(robot,[path(i,:) theta(i)]);
    end
    hold off

    more_qf=input('enter 0 to stop, 1 to continue: ');

end

prmpath=path';
waypoints=vecnorm(diff(prmpath')');
waypointslength=[0 cumsum(waypoints)];
pathlength=sum(waypoints);
pathvel=1; % m/s
ts=.5;
lambda=(0:pathvel*ts:pathlength);n_l=length(lambda);
despath=zeros(2,n_l);
despath(1,:)=interp1(waypointslength,prmpath(1,:),lambda);
despath(2,:)=interp1(waypointslength,prmpath(2,:),lambda);
figure(50);plot(prmpath(1,:),prmpath(2,:),'x-',despath(1,:),despath(2,:));
view(-90,90);axis([-1 11 -1 11]);
save prmpath despath


