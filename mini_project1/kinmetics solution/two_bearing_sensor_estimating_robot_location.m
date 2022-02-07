%
% illustration of using two bearing sensors 
% (how do you find center of the circle?)
% 
%clear all;close all;
pB1=[1;2];
pB2=[4;1];
l=norm(pB1-pB2);
e12=(pB1-pB2)/l;

% assume bearing sensor outputs phi = phi_1 - phi_2
%phi=pi/3;
phi=rand*pi;
% if phi is negative then switch pB1 and PB2 then phi is positive

% try different psi angles (angle from pB2-pB1 to pB2-pR)
% from 0 to pi-phi
N=100;
psi=(pi/N:pi/N:pi-phi-pi/N);
N=length(psi);

for i=1:N
    % k2 is the unit vector from pB2 to pB1 rotated by psi
    % k1 is the unit vector from pB1 to pB2 rotated by -(pi-phi-psi)
    k2=rot2(psi(i))*e12;k1=rot2(-(pi-phi-psi(i)))*(-e12);
    % solve for the line intersections 
    alpha(:,i)=inv([-k1 k2])*(pB1-pB2);
    % the intersection is the robot position
    pR(:,i)=pB1+alpha(1,i)*k1;
    % plot to show intersection of the two lines
%     a=(0:.1:5);
%     figure(20);
%     plot(pB1(1)+k1(1)*a,pB1(2)+k1(2)*a,pB2(1)+k2(1)*a,...
%         pB2(2)+k2(2)*a,pB1(1),pB1(2),'x',pB2(1),pB2(2),'o',...
%         pR(1,i),pR(2,i),'s','LineWidth',2)
    %pause(1);
end

% show all the intersections
figure(1);plot(pR(1,:),pR(2,:),'x',[pB1(1),pB2(1)],[pB1(2) pB2(2)],'o',...
    'linewidth',2);
title(['\phi=',num2str(phi*180/pi),' deg']);
axis([-1 5 -3 3]);axis('square');
hold on

%
for i=1:N
    plot([pR(1,i) pB1(1)],[pR(2,i) pB1(2)],'k-',...
        [pR(1,i) pB2(1)],[pR(2,i) pB2(2)],'r-')
end

% check indeed phi is what it started out to be

z1=(pB1-pR)./vecnorm(pB1-pR);
z2=(pB2-pR)./vecnorm(pB2-pR);
for i=1:N;phicheck(i)=subprob1([0;0;1],[z2(:,i);1],[z1(:,i);1]);end
disp('difference between solution phi vs. actual phi')
display(norm(phicheck-phi));

% find the circle 
[r,pC,k]=threepointscircle([pB1 pB2 pR(:,1);zeros(1,3)]);

M=200;th=(0:2*pi/M:2*pi*(1-1/M));
pCircle=pC(1:2)+[r*cos(th);r*sin(th)];
figure(1);plot(pCircle(1,:),pCircle(2,:),':','LineWidth',3);

phi1=(rand-.5)*2*pi;phi2=phi1-phi;
[pC1,r1,q1,qrange]=two_bearings([pB1 pB2],[phi1;phi2]);

M=30;th=(q1:qrange/M:q1+qrange);
figure(1);plot(pC1(1)+r1*cos(th),pC1(2)+r1*sin(th),'^','LineWidth',3);

hold off
