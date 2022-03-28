%
% 2D formation control example based on the relative pose
% 
clear all;close all;
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
        zstar1 = [-.2;-.2]*3;
        zstar2 = [.2;-.2]*3;
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
        zc = 2;
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
% # of time steps
N = 400;
% allocate space
q = zeros(n*p,N+1);
u = zeros(n*p,N);
z = zeros(n*l,N+1);
% initial condition
q(:,1) = rand(n*p,1);
z(:,1) = kron(D',eye(n,n))*q(:,1);
% time step
ts = 1;
% time vector
t = (0:N)*ts;
% feedback gain
% for collision avoidance with each other
Kp_repel = 0.15*eye(n*p,n*p);
% for keeping the formation together
Kp_form = 0.5*eye(n*p,n*p);
% define repellent potential function 
% (only active for distance within qrange) 
qrange=.05;
psifun = @(x,dx) (x<dx).*(1./x-1./dx);
% desired motion for the formation 
% rc = center of motion, w = angular frequency
% vmag = magnitude for the formation velocity
% th = initial target for the lead agent
% kp1 = gain for the lead agent
%rc = 4;w=.05;vmag=.1;th=pi/4;kp1=0.3;
%rc = 4;w=.05;vmag=0;th=pi/4;kp1=0;
rc = 2;w=.05;vmag=0;th=pi/4;kp1=0.3;
%rc = 4;w=.05;vmag=2;th=pi/4;kp1=0;
% ****************
% control loop
% ****************
for k=1:N    
    % set up desired formation motion target at time t(k)
    vd = vmag * w *[-sin(w*t(k));cos(w*t(k))];
    vd = reshape(vd*ones(1,p),n*p,1);
    qd1 = rc*[cos(w*t(k)+th);sin(w*t(k)+th)];
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
    % add in lead agent motion if specified (depending on kp1)
    %u(1:n,k) = u(1:n,k)-kp1*(q(1:n,k)-qd1);
    u(:,k) = u(:,k) + kron(ones(p,1),-kp1*(q(1:n,k)-qd1));
    % update agent position
    q(:,k+1) = q(:,k)+ts*u(:,k);
    % update formation difference variable
    z(:,k+1) = kron(D',eye(n,n))*q(:,k+1);
end

% plot motion of agents in the plane and save as a movie
figure(2);
set(gcf,'position',[500,400,400,400]);
for k=1:N+1
    qk=q(:,k);
    qkflat=reshape(qk,n,p);
    plot(qkflat(1,1),qkflat(2,1),'r^',...
        qkflat(1,2:p),qkflat(2,2:p),'b^','linewidth',3);
    axis([-10,10,-10,10]);axis('square');
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
