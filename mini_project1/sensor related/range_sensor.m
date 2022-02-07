%
% range_sensor.m
%
% input:    pL = 3xN landmark locations
%           pR = 3x1 robot location
%           sigma = zero mean Gaussian noise standard deviation
%

function y=range_sensor(pL,pR,sigma)

y=vecnorm(pL-pR)'+randn(size(sigma)).*sigma;

end