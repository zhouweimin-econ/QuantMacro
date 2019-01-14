function z = scaleDown(x,xmin,xmax)

% This function maps x in [xmin,xmax] into z in [-1,1]; inverse of scale_up
% 
% Inputs
%   (1) x: array with elements in [xmin,xmax]; if outside range, will force
%   (2) xmin: scalar lower bound of domain 
%   (3) xmax: scalar upper bound of domain
%
% Outputs
%   (1) z: array with elements in [-1,1]
% 
% Thomas Winberry, January 19, 2016

z = min(max(2 * ((x - xmin) / (xmax - xmin)) - 1,-1 * ones(size(x))),ones(size(x)));