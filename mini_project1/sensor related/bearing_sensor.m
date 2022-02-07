%
% bearing_sensor.m
%
% input:    pB = 3xN landmark locations
%           k_heading = 3x1 unit vector of robot's heading direction
%           k_vertical = 3x1 unit vector of sensor's normal axis
%           (usually perpendicular to k_heading)
%           sigma = zero mean Gaussian noise standard deviation
%

function y=bearing_sensor(pB,k_heading,k_vertical,sigma)

k_heading=k_heading/norm(k_heading);

for i=1:size(pB,2)

    y(i)=subprob1(k_vertical,k_heading,pB(:,i)/norm(pB(:,i))) + ...
        randn*sigma(i);

end

end

function khat = hat(k)
  
  khat=[0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];

end