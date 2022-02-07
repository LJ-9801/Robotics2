%
% use keyboard to drive a wheeled mobile robot
% 
% to run: just type in 
% 
% >> autodrive_scan
% 
% use the arrow key in the figure window to drive the robot 
% (up=forward, down=backware, left=turn left, right=turn right)
% press q to quit
%
% this program calls 
%
% roomspec: set up the room and all collision objects in colobj structure
% roomshow: show all the collision objects
% robotspec: set up the robot collision body
% robotshow: show robot in the same plot
% robotfcn: function that gets called when a key is pressed in the figure
% output: generates all the sensor reading based on the current robot pose
% pose_est: estimating the pose of the robot
% colcheck: collision check between robot and all objects in room
%

% input u is assume available

% define room
% roomspec;
% % show room
% fignum=1;
% h_fig=roomshow(colobj,fignum);
% axis('square');
% % define robot
% rL=.4;rW=.2;rz=.1;
% robot=robotspec([rL;rW;2*rz]);
% % show robot
% q0=[.2;2;0];
% q(:,1)=q0;
% robotshow(robot,q0);
% %
% ez=[0;0;1];
% **** range sensors ****
% UWB (local GPS)
% zW=colobj.obj{1}.Z;
% lW=colobj.obj{1}.X;
% zW=0;
% pL(:,1)=[0;0;zW];pL(:,2)=[lW;0;zW];pL(:,3)=[0;lW;zW];pL(:,4)=[lW;lW;zW];
% % **** bearing sensors ****
% pB(:,1)=colobj.obj{14}.Pose(1:3,4);
% pB(:,2)=colobj.obj{10}.Pose(1:3,4);
% %
% UWB (local GPS)
% zW=colobj.obj{1}.Z;
% lW=colobj.obj{1}.X;
% zW=0;
% pL(:,1)=[0;0;zW];pL(:,2)=[lW;0;zW];pL(:,3)=[0;lW;zW];pL(:,4)=[lW;lW;zW];
% **** bearing sensors ****
%pB(:,1)=colobj.obj{14}.Pose(1:3,4);
%pB(:,2)=colobj.obj{10}.Pose(1:3,4);
%pB(3,:)=[0 0];
%N_scan=8;
%N_range=size(pL,2);N_bearing=size(pB,2);N_odo=2;
%
%ns=N_range+N_bearing+N_odo+N_scan; % total # of sensors
%wcov=[0.05;0.05];vcov=.15*ones(ns,1); % noise covariance
% steering command and sampling period
%v=.1;w=.1;ts=1;
% initial sensor reading
%y(:,1)=output(q(:,1),pL,pB,N_scan,[0;0],vcov,robot,colobj);
% initial state estimate
% qhat(:,1)=pose_est(y(:,1),pL,pB,N_scan,wcov,vcov);
% initialization for the EKF
qhat_EKF(:,1)=[0;0;0];P{1}=.1*eye(3,3);S{1}=.1*eye(3,3);

N_step=size(u,2);
%N_step=50;
% scan output range
% # of scan lines in lidar
ns=(N_range+N_odo+N_bearing+1:N_range+N_odo+N_bearing+N_scan);
% figure #
fignum=10;
% show one scan in robot frame
figure(fignum);
% no initial knowledge
qhat(:,1)=[0;0;0];
% first scan
[xloctemp,yloctemp]=scanloc(y(ns,1),zeros(3,1));
xloc_k(:,1)=xloctemp+qhat(1,1);
yloc_k(:,1)=yloctemp+qhat(2,1);
plot(xloctemp,yloctemp,'x');
hold on;
view(-90,90);axis([-1 11 -1 11 0 4])
axis([-1 11 -3 9 0 4])
% up to max = 3 segments
numseg=2;
% minimum segment length
minL=3;
% alignment threshold 
threshold=.2;
% extract lines
%[Pts,SegLength,Idx,unitvec,mbline,xline,yline]=...
%    linesegment([xloctemp yloctemp],threshold,numseg,minL,30);
%

N_step=size(u,2);
for k=1:N_step
    % generate output
    y(:,k+1)=output(q(:,k+1),pL,pB,N_scan,utrue(:,k),vcov,robot,colobj);
    % kinematics 
    [qhat(:,k+1),uu]=wmr(qhat(:,k),utrue(:,k),ts,zeros(size(wcov)));
    % one scan in robot frame
    [xloctemp,yloctemp]=scanloc(y(ns,k),zeros(3,1));
    % estimate robot state (wrt initial robot pose) using the robot
    % kinematics 
    % EKF estimate
%     [qhat_EKF(:,k+1),P{k+1},S{k+1}]=... 
%         pose_est_kalman(qhat_EKF(:,k),u(:,k),y(:,k),ts,...
%         pL,pB,N_scan,wcov,vcov,P{k},S{k},robot,colobj);
    % switch to the estimate frame
    scanloc_k=(rot2(qhat(3,k+1))*[xloctemp yloctemp]'+qhat(1:2,k+1))';
    % extract lines
    [Pts,SegLength,Idx,unitvec,mbline,xline,yline]=...
        linesegment(scanloc_k,threshold,numseg,minL,30);
    %[Pts,Idx,lineparam,err,xline,yline]=extractline(scanloc_k);
    % 
    figure(fignum);plot(scanloc_k(:,1),scanloc_k(:,2),'.');
    figure(fignum);plot(xline{1},yline{1},'b-');
    disp([xline{1};yline{1}]);
    xloc_k(:,k+1)=scanloc_k(:,1);
    yloc_k(:,k+1)=scanloc_k(:,2);
end

hold off;
