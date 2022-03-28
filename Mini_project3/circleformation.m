clear all;close all;

% arbitrary # of agents
p = 10;
%
% selection 
graphsel=5; % default
graphsel=input(['enter 1 for 3 nodes in V configuration \n',...
    '2 for 3 nodes fully connected \n',...
    '3 for a random number of nodes in a ring\n',...
    '4 for a random number of nodes fully connected\n',...
    '5 for a random number of nodes in random configuration:\n']);
graphsel=min(max(graphsel,1),5);
switch graphsel
    case 1
        % 3 agent V-shaped
        D = [-1 -1;1 0;0 1];B=diag([2 1 1]);
        G = graph(B-D*D');
    case 2
        % 3 agent fully connected
        D = [-1 -1 0 ; 1 0 -1; 0 1 1];B=diag([2 2 2]);
        G = graph(B-D*D');
    case {3,4,5}
        p=input('enter # of agents (or 0 for random): ');
        if p<=0;p=fix((rand+.5)*30);end
        % p agents
        switch graphsel
            case 3
                % p agent in a ring
                n1=(1:p);
                n2=[(2:p) 1];
            case 4
                % p agents fully connected
                n1=[];n2=[];
                for i=1:p
                    n1=[n1 i*ones(1,p-i)];
                    n2=[n2 (i+1:p)];
                end
            case 5
                % random p agent connection
                nodepair=[randperm(2*p);randperm(2*p)];
                nodepair=mod(nodepair,p)+1;
                ind=(nodepair(1,:)==nodepair(2,:));
                nodepair(:,ind)=[];
                ns=sortrows(sort(nodepair)')';
                vecnorm(diff(ns')');
                ind1=find(vecnorm(diff(ns')')>0);
                ind1=[ind1(1) ind1+1];
                ns1=ns(:,ind1);
                n1=ns1(1,:);n2=ns1(2,:);
        end
        % generate graph using node pairs
        G=graph(n1,n2);
        % generate the non-sparse incidence matrix from the graph
        D=full(incidence(G));
end
% show the graph connection of agents
figure(77);plot(G);
set(gcf,'position',[100,400,500,500]);
% # of links
l = size(D,2);
% # of states for each agent
n = 2;
% center for formation
c = [0;0];
% # of time steps
N = 100;
% # of agents
p = size(D,1);
% allocate space
q = zeros(n*p,N+1);
u = zeros(n*p,N);
z = zeros(n*l,N+1);
% initial condition
q(:,1) = rand(n*p,1);
z(:,1) = kron(D',eye(n,n))*q(:,1);
% distance from center
dc = .2;
% distance from each other
%dz = [2;2;2];
dz = dc*2*pi/p*ones(l,1);
dz = 2*ones(l,1);
% time step
ts = 1;
% time vector
t = (0:N)*ts;
% feedback gain
Kp = .1*eye(n*p,n*p);
Kp1 = 0.1*eye(n*p,n*p);
% define the potential function
%s = @(x,d) (1/d - 1/x);
s = @(x,d) log(x/d);

% control loop
for k=1:N
    % link difference vectors z = D^T * q 
    zk = kron(D',eye(n,n))*q(:,k);
    % for each link compute the potential gradient
    for j=1:l
        zkj = zk((j-1)*n+1:j*n);        
        psi((j-1)*n+1:j*n,1) = s(norm(zkj),dz(j))*zkj/norm(zkj);
    end
    % potential gradient to the origin
    for i=1:p
        qki = q((i-1)*n+1:i*n,k);
        phi((i-1)*n+1:i*n,1) = s(norm(qki-c),dc)*(qki-c)/norm(qki-c);
    end
    % controller combines the two
    u(:,k) = -Kp*kron(D,eye(n,n))*psi - Kp1 * phi;
    q(:,k+1) = q(:,k)+ts*u(:,k);
    z(:,k+1) = kron(D',eye(n,n))*q(:,k+1);
end



% figure(1);hold on;
% for i=1:p
%     plot(q((i-1)*n+1,:),q(i*n,:),'k.','linewidth',2);
%     %plot(q((i-1)*n+1,1),q(i*n,1),'bo','linewidth',3);
%     plot(q((i-1)*n+1,N+1),q(i*n,N+1),'rx','linewidth',3);
%     if i<p
%         plot([q((i-1)*n+1,N+1),q(i*n+1,N+1)],...
%             [q(i*n,N+1),q((i+1)*n,N+1)],'c-','linewidth',3);
%     end
%     plot(c(1),c(2),'ys','linewidth',5);
%     axis('square');grid;
% end
zf = kron(D',eye(n,n))*q(:,N+1);
disp('final z distance');
disp(vecnorm(reshape(zf,n,l)))
% qf= reshape(q(:,N+1),n,p);
% disp('final distance from c');
% disp(vecnorm(qf-c));
% hold off

figure(2);
for k=1:N+1
    qk=q(:,k);
    qkflat=reshape(qk,n,p);
    plot(qkflat(1,:),qkflat(2,:),'^','linewidth',3);
    axis([-4,4,-4,4]);axis('square');
    M(k)=getframe;
end
movie(M);

figure(10);hold on
for i=1:l
    zvec{i}(:,:)=z((i-1)*n+1:i*n,:);
    plot(t,vecnorm(zvec{i})-dz(i));
end
hold off

