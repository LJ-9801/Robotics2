%
% threeballs.m
% 
% find intersections of three balls with centers at p1,p2,p3 
%   and radii r1, r2, r3
% returns up to two solutions 
% 
function z=threeballs(p1,p2,p3,r1,r2,r3)

% find the unit vector k from ball 1 center to ball 2 center
k=(p2-p1)/norm(p2-p1);

% find a vector at the intersection of ball 1 and ball 2
% call this vector r1rot
% 1. first find the unit vector perpendicular to both k and p1 (kperp)
kperp=randn(3,1);kperp=kperp-kperp'*k*k;kperp=kperp/norm(kperp);
% 2. find the intersection of the cone rotating r1*k about kperp that
% intersects ball 2.  this means solving for q1 in 
% ||R(kperp,q1)*r1*k + p1 - p2|| = r2
%
q1=subprob3(kperp,r1*k,p2-p1,r2);
r1rot=rot(kperp,q1(1))*k*r1;
% now look for the intersection of cone rotating r1rot about k and ball 3
q=subprob3(k,r1rot,p3-p1,r3);

% there are two solutions
z(:,1)=rot(k,q(1))*r1rot+p1;
z(:,2)=rot(k,q(2))*r1rot+p1;

end
