%
% 2D formation control example based on the relative pose in a room
% 
clear all;close all;

% *****
% define room and robot
% *****
% define room
roomspec;
% define robot
rL=.4;rW=.2;rz=.1;
robot=collisionCylinder(rW,rz);
% *****
% import desired path (assuming same sampling period ts)
% variable name of path: despath
% *****

while ~exist('despath');
    filename=input('enter desired path file name with variable despath: ','s');
    eval(['load ',filename]);
    if ~exist('despath');disp('No variable named *despath*');end
end

% # of time steps
N = length(despath);

% *****
% formation control part
% *****
% repeatable runs
rng(100);
% ask user for # of agents
p = 13;
p = input('enter # of agents: ');
p=p+~mod(p,2);
%
% ask user for the formation type
%
formation_type = input('enter 1 for V, 2 for circle formation: ');
switch formation_type
    % 1 = V formation with 1=lead node 2,..,(p-1)/2 for right branch,
    % (p-1)/2+1,..,p for left branch
    case {1}
        % set up the incidence matrix for V graph
        D=zeros(p,p-1);
        D(1,1)=-1;D(1,1+(p-1)/2)=-1;
        D((p-1)/2+1,(p-1)/2)=1;D(p,(p-1))=1;
        for i=2:(p-1)/2
            D(i,i-1)=1;D(i,i)=-1;
            D(i+(p-1)/2,(p-1)/2+i-1)=1;D(i+(p-1)/2,(p-1)/2+i)=-1;
        end
        % set up the target difference variables
        zstar1 = [-.3;-.3];
        zstar2 = [.3;-.3];
        zstar = [kron(ones((p-1)/2,1),zstar1);kron(ones((p-1)/2,1),zstar2)];
    % 2 = circle formation 
    otherwise
        % set up incidence matrix for a ring graph
        D=zeros(p,p);
        D(1,1)=-1;D(1,p)=1;
        D(p,p-1)=1;D(p,p)=-1;
        for i=2:p-1
            D(i,i-1)=1;D(i,i)=-1;
        end
        % set up target difference variables 
        zc = .8;
        zstar = zeros(2*p,1);
        zstar(1:2)=zc*[cos(pi/p);-sin(pi/p)];
        for i=2:p
            zstar((i-1)*2+1:2*i)=...
                zc*[cos((i-1)*2*pi/p+pi/p);-sin((i-1)*2*pi/p+pi/p)];
        end
end
% show graph 
G=graph(-D*D'+diag(diag(D*D')));
figure(20);plot(G);
set(gcf,'position',[100,400,400,400]);
%
% for collision avoidance between agents, use fully connected graph
%
Dfull=(incidence(graph(ones(p,p)-diag(ones(p,1)))));
% *************************************
% # of links for the formation graph
l = size(D,2);
% # of states for each agent
n = 2;
% allocate space
q = zeros(n*p,N+1);
u = zeros(n*p,N);
z = zeros(n*l,N+1);
% initial condition
q(:,1) = rand(n*p,1);

q0=[.4;2;0];
q0=[rW+.1;2;0];
q0flat=zeros(n,p);
for i=1:p;q0flat(:,i)=[rW+.1;2.2*rW*i];end
q(:,1)=reshape(q0flat,n*p,1);
z(:,1) = kron(D',eye(n,n))*q(:,1);
% time step
ts = 1;
% time vector
t = (0:N)*ts;
% feedback gain
% for collision avoidance with each other
Kp_repel = 0.15*eye(n*p,n*p);
% for keeping the formation together
Kp_form = 0.2*eye(n*p,n*p);
% define repellent potential function 
% (only active for distance within qrange) 
qrange=.1;
psifun = @(x,dx) (x<dx).*(1./x-1./dx);
% lead agent feedback gain
Kp1 = .2;
% obstacle avoidance gain
k_obs = .0008;
% *****
% Show robot and room
% *****
figure(1);hold on;
roomshow(colobj,1);axis('square');
view(-90,90);axis([-1 11 -1 11 0 4])
for i=1:p;robotshow(robot,[q0flat(:,i);0],[.2,.1,.5]);end
hold off

% ****************
% control loop
% ****************
for k=1:N    
    % set up desired formation motion target at time t(k)
    if k<N
        vd = kron(ones(p,1),despath(:,k+1)-despath(:,k));
    else
        vd = zeros(size(vd));
    end
    % ****
    % formation control part
    % ****
    % link difference vectors z = D^T * q 
    zk = kron(D',eye(n,n))*q(:,k);    
    % gradient to the formation
    phi = kron(D,eye(n))*(zk-zstar);
    % ****
    % repellent function to avoid inter-agent collision
    % ****
    zkfull=kron(Dfull',eye(n,n))*q(:,k);
    zkflat=reshape(zkfull,n,size(Dfull,2));
    distzk=vecnorm(zkflat);
    % if within qrange, push each other away along the distance vector
    psiflat = zkflat.*psifun(distzk,qrange)./distzk;
    psi = kron(Dfull,eye(n,n))*reshape(psiflat,n*size(Dfull,2),1);
    % ****
    % controller combines formation and repellent parts and add in desired
    % formation motion velocity
    % ****
    u(:,k) = -Kp_repel*psi - Kp_form * phi + vd;
    % lead agent knows the path
    %u(1:2,k) = u(1:2,k) - Kp1 * (q(1:2,k)-despath(:,k));
    u(:,k) = u(:,k) + kron(ones(p,1),- Kp1*(q(1:2,k)-despath(:,k)));    
    % *****
    % obstacle aovidance 
    % *****
    % check current step
    qkflat=reshape(q(:,k),n,p);
    for j=1:p
        [isInt,dist,wp]=colcheck(robot,[qkflat(:,j);0],colobj);
        [rho,ind]=min(dist);
        dwp=wp{ind}(1:2,2)-wp{ind}(1:2,1);
        u(n*(j-1)+1:n*j,k)=u(n*(j-1)+1:n*j,k) - ...
            k_obs*psifun(rho,qrange)*dwp/rho;
        if(max(isnan(dist))>0);
            fprintf('collision, k=%d, i=%d\n',k,i);
        end
    end
    % update agent position
    q(:,k+1) = q(:,k)+ts*u(:,k);
    %
    % check next step
    %
    qkp1flat=reshape(q(:,k+1),n,p);
    for j=1:p
        [isInt,dist,wp]=colcheck(robot,[qkp1flat(:,j);0],colobj);
        while (max(isnan(dist)>0))
            qkp1flat(:,j)=qkflat(:,j)+randn(2,1)*.01;
            [isInt,dist,wp]=colcheck(robot,[qkp1flat(:,j);0],colobj);
        end
    end
    %
    q(:,k+1)=reshape(qkp1flat,n*p,1);
    % *****
    % update formation difference variable
    z(:,k+1) = kron(D',eye(n,n))*q(:,k+1);    
%     figure(1);hold on
%     roomshow(colobj);axis('square');
%     view(-90,90);axis([-1 11 -1 11 0 4])
%     for i=1:p;robotshow(robot,[qkflat(:,i);0],[.5,.1,.5]);end
%     hold off;
end

% plot motion of agents in the plane and save as a movie
figure(2);
set(gcf,'position',[500,400,400,400]);
for k=1:N+1
    qk=q(:,k);
    qkflat=reshape(qk,n,p);
    figure(2);hold on
    roomshow(colobj,1);axis('square');
    view(-90,90);axis([-1 11 -1 11 0 4])
    for i=1:p;robotshow(robot,[qkflat(:,i);0],[.8,.1,.1]);end
    hold off;
    M(k)=getframe;
end
%movie(M);

% check the distance between agents and show minimum distance
zfull=kron(Dfull',eye(2))*q;
for k=1:N+1;
    minznorm(k)=min(vecnorm(reshape(zfull(:,k),n,size(Dfull,2))));
end
figure(3);plot(t,minznorm,'LineWidth',2);grid;
set(gcf,'position',[100,1,400,400]);
xlabel('time (s)');ylabel('minimum agent distance (m)');
disp('minimum agent distance')
disp(min(minznorm));

% plot the difference variable (should converge to zero) to see how
% formation holds during motion
figure(10);hold on
set(gcf,'position',[500,1,400,400]);
dz=z-zstar;
for i=1:l
    plot(t,vecnorm(dz((i-1)*n+1:i*n,:)),'linewidth',2);
    xlabel('time (s)');ylabel('formation difference error');
end
hold off
% show final formation difference variable
disp('final delta z');
disp((reshape(dz(:,end),n,l)))
