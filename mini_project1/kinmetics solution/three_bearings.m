%
% three_bearing.m
%
% input:
%   pB = bearing locations (in plane) (2x3 vector, column i = pB_i)
%   phi =  bearing angle measurements (3x1 vector)
%
% output:
%   pR = solution of robot location 
%   qR = robot orientation angle
%

function [pR,qR]=three_bearings(pB,phi)

% use two bearing sensors to find the soution circles
[pCA,rA,qA1,qArange]=two_bearings(pB(1:2,1:2),phi(1:2));
[pCB,rB,qB1,qBrange]=two_bearings(pB([1,3],[1,3]),phi([1,3]));

% find intersections of the two circles
ex=[1;0;0];ez=[0;0;1];
qA=subprob3(ez,rA*ex,[pCB-pCA;0],rB);
%qB=subprob3(ez,rB*ex,[pCA-pCB;0],rA);
kAB=(pCB-pCA)/norm(pCB-pCA);
qAB=subprob1(ez,ex,[kAB;0]);
qB(1)=(pi-asin((rA/rB)*sin(qA(1)-qAB)))+qAB;
qB(2)=(pi-asin((rA/rB)*sin(qA(2)-qAB)))+qAB;
qA=anglerange(qA);
qB=anglerange(qB);
qA=(qA<0).*(qA+2*pi)+(qA>0).*qA;
qB=(qB<0).*(qB+2*pi)+(qB>0).*qB;

% use only the one that is within both allowable ranges

qsol=[];
if ((qA(1)>qA1+100*eps)&(qA(1)-qA1<=qArange))&...
        ((qB(1)>qB1+100*eps)&(qB(1)-qB1<=qBrange))
    qsol=qA(1);
end

if ((qA(2)>qA1+100*eps)&(qA(2)-qA1<=qArange))&...
        ((qB(2)>qB1)&(qB(2)-qB1<=qBrange+eps))
    qsol=[qsol;qA(2)];
end

%disp([qsol',qA'])
qsol=qA;

for i=1:length(qsol)
    % find pR
    pR(:,i)=pCA+rot2(qsol(i))*[rA;0];

    % then qR
    qR(i)=subprob1(ez,ex,[rot2(-phi(1))*(pB(1:2,1)-pR(:,i))...
        /norm(pB(1:2,1)-pR(:,i));0]);
end

end