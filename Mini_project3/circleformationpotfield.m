%
% circleformationpotfield.m
%
% using potential field to keep agents at set distance apart around a
% circle
% 
% Select option 4 to see how the circle formation using distance metric 
%
clear all;close all;
%

% selection 
graphsel=4; % default
graphsel=input(['enter 1 for 3 nodes in V configuration \n',...
    '2 for 3 nodes fully connected \n',...
    '3 for a random number of nodes in a ring\n',...
    '4 for a random number of nodes fully connected\n',...
    '5 for a random number of nodes in random configuration:\n']);
graphsel=min(max(graphsel,1),5);

% 
switch graphsel
    case 1
        % V-shaped
        D = [-1 -1;1 0;0 1];B=diag([2 1 1]);
        G = graph(B-D*D');
    case 2
        % fully connected
        D = [-1 -1 0 ; 1 0 -1; 0 1 1];B=diag([2 2 2]);
        G = graph(B-D*D');
    case {3,4,5}
        p=input('enter # of agents (or 0 for random): ');
        if p<=0;p=fix((rand+.5)*30);end
        switch graphsel
            case 3
                % ring
                n1=(1:p-1);
                n2=[(2:p)];
            case 4
                % fully connected
                n1=[];n2=[];
                for i=1:p
                    n1=[n1 i*ones(1,p-i)];
                    n2=[n2 (i+1:p)];
                end
            case 5
                % random graph (may not be connected)
                % the code below captues all the nodes and remove
                % duplicates
                np=fix(1.1*p);
                nodepair=[randperm(np);randperm(np)];
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
        G=graph(n1,n2);
        D=full(incidence(G));
end
figure(20);plot(G);
set(gcf,'position',[100,400,500,500]);
% # of links
l = size(D,2);
% # of states for each agent
n = 2;
% # of time steps
N = 500;
% # of agents
p = size(D,1);
% allocate space
q = zeros(n*p,N+1);
u = zeros(n*p,N);
% initial condition
q(:,1) = rand(n*p,1);
% circular path for the desired formation motion
lambda = (0:N)/N;
th = lambda*2*pi;
rc = 2.5;
ctraj = rc*[cos(th);sin(th)];
% distance from center
dc = 2;
% distance from each other
% spread out around the circle
znormstar = dc*2*pi/p*ones(l,1);
% time step
ts = 1;
% time vector
t = (0:N)*ts;
% feedback gain
K_repel = .2*eye(n*p,n*p); % repelling gain
K_track = 2.8*eye(n*p,n*p); % tracking gain
% define the potential function
% s = @(x,d) (1/d - 1/x);
s = @(x,d) log(x/d);
%
qrange=1.2*znormstar(1);
psifun = @(x,dx) (x<dx).*(1./x-1./dx);

for k=1:N
    qk=q(:,k);
    qkflat=reshape(qk,n,p);
%    
    zk=kron(D',eye(n,n))*qk;
    zkflat=reshape(zk,n,l);
    distzk=vecnorm(zkflat);
    psiflat = zkflat.*psifun(distzk,qrange)./distzk;
    psi = reshape(psiflat,n*l,1);
%    
    for i=1:p
        qki=qkflat(:,i);
        c=ctraj(:,k);
        phi((i-1)*n+1:i*n,1) = s(norm(qki-c),dc)*(qki-c)/norm(qki-c);
    end
    u(:,k) = K_repel*kron(D,eye(n,n))*psi - K_track*phi;
    q(:,k+1) = q(:,k)+ts*u(:,k);
end

zf=reshape(kron(D',eye(n,n))*q(:,N+1),n,l);
disp('final separation');
disp(vecnorm(zf));
figure(15);plot(vecnorm(zf),'o','linewidth',2)

figure(1);hold on;
for i=1:p
    plot(q((i-1)*n+1,:),q(i*n,:),'k.-','linewidth',2);
    plot(q((i-1)*n+1,1),q(i*n,1),'bo','linewidth',3);
    plot(q((i-1)*n+1,N+1),q(i*n,N+1),'rx','linewidth',3);
    if i<p
        plot([q((i-1)*n+1,N+1),q(i*n+1,N+1)],...
            [q(i*n,N+1),q((i+1)*n,N+1)],'c-','linewidth',3);
    end
    plot(c(1),c(2),'ms','linewidth',5);
end
qf= reshape(q(:,N+1),n,p);
disp('final distance from c');
disp(vecnorm(qf-c));
hold off
axis([-2 4 -2 4]);axis('square');
figure(2);
for k=1:N+1
    qk=q(:,k);
    qkflat=reshape(qk,n,p);
    plot(qkflat(1,:),qkflat(2,:),'^',...
        ctraj(1,k),ctraj(2,k),'ms','linewidth',3);
    axis([-4,4,-4,4]);axis('square');
    M(k)=getframe;
end
%movie(M);

