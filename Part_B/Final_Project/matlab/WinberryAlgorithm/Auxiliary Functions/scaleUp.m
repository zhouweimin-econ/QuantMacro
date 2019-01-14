function x = scaleUp(z,xmin,xmax)

% This function maps z in [-1,1] into x in [xmin,xmax]; inverse of scale_down
% 
% Inputs
%   (1) z: array with elements in [-1,1]; if outside range, will force
%   (2) xmin: scalar lower bound of range
%   (3) xmax: scalar upper bound of range
%
% Outputs
%   (1) x: array with elements in [xmin,xmax]
% 
% Thomas Winberry, January 19, 2016

x = min(max((.5 * (z + 1) * (xmax - xmin)) + xmin,xmin * ones(size(z))),...
    xmax * ones(size(z)));
