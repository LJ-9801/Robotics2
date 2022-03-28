%
% example of phase synchronization using phase difference (as a graph)
% feedback 
%
% choose a case to visualize
%   1 = 3 nodes in V graph
%   2 = 3 nodes in fully connected graph (triangle)
%   3 = random # of nodes in a ring
%   4 = random # of nodes fully connected
%   5 = random # of nodes in random connectivitiy
%

clear all;close all;
%
upperright=[800 100 400 400];
leftside=[100 400 400 400];

% selection 
graphsel=5; % default
graphsel=input(['enter 1 for 3 nodes in V configuration \n',...
    '2 for 3 nodes fully connected \n',...
    '3 for a random number of nodes in a ring\n',...
    '4 for a random number of nodes fully connected\n',...
    '5 for a random number of nodes in random configuration:\n']);
graphsel=min(max(graphsel,1),5);
p=fix((rand+.5)*30);

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
set(gcf,'position',upperright);

% # of links
l=size(D,2);
% # of nodes
n=size(D,1);
% feedback gain 
Kp=diag((l:-1:1))*.2;
% 
omega0=0.8; % Hz
% full feedback u = omega0 - D*Kp*D^T

% use weighted Laplacian feedback
% add a state to acount for omega_0
A=[-D*Kp*D' ones(n,1);zeros(1,n+1)];
tmax=20;
t=(0:.1:tmax);

% random initial phase
q0=(rand(n,1)-.5)*2*pi;

% solve for the closed loop system
for i=1:length(t)
   x = expm(A*t(i))*[q0;omega0];
   q(:,i)=x(1:n);
   om(i)=x(n+1);
end

% plot the phases
figure(1); h1=plot(t,q,'linewidth',2);
set(gcf,'Position',leftside);
title('phase angles');xlabel('time (sec)');ylabel('\theta');

% plot as phase of sinusoids (shifted vertically for visualization)
figure(3);hold on;
for i=1:n
    plot(t,sin(q(i,:))+2*(i-1),'linewidth',2);
end
hold off 

% show phase as rotating arrow

for k=1:length(t)
    figure(4);
    arrow3(zeros(n,2),[cos(q(:,k)) sin(q(:,k))]);
    axis([-1 1 -1 1]);axis('square');
    pause(.01);
end
