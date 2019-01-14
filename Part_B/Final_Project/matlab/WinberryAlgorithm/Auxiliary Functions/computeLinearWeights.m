function [vIndicesBelow,vIndicesAbove,vWeightBelow,vWeightAbove] = ...
	computeLinearWeights(vGrid,vValues);

% Given a grid and values off the grid, finds (1) points above and below on the grid and (2) distance
% from these points.  Useful for linear interpolation and constructing histogram transition.
%
% Inputs
%	(1) vGrid: grid of points find distances between (column vector)
%	(2) vValues: vector of values off the grid, to compute weights to closest grid points (column vector)
%	
% Outputs
%	(1) vIndicesBelow: for each point in vValues, index of gridpoint immediately below
%	(2) vIndicesAbove: for each point in vValues, index of gridpoint immediately above
%	(3) vWeightBelow: for each point in vValues, computes distance between value and grid point
%		above, and then expresses as fraction of total size of bracketing interval (used
%		in constructing linear interpolation and histogram transition).  So if value very close to
%		lower gridpoint, gets more weight.
%	(3) vWeightAbove: for each point in vValues, computes distance between value and grid point
%		below, and then expresses as fraction of total size of bracketing interval (used
%		in constructing linear interpolation and histogram transition)
%
% Thomas Winberry, July 26th, 2016

% Find nearest point to vValues on the asset grid
vIndices = knnsearch(vGrid,vValues);		% find index of nearest neighbor on grid
vGridIndices = vGrid(vIndices);	% find value of nearest neighbor on grid

% Find indices above and below choice
vIndicesBelow = vIndices;  vIndicesAbove = vIndices; [nGrid,~] = size(vGrid);

vIndicesBelow(vGridIndices > vValues) = vIndicesBelow(vGridIndices > vValues) - 1;
vIndicesBelow(vIndicesBelow <= 1) = 1;
vGridIndicesBelow = vGrid(vIndicesBelow);

vIndicesAbove(vGridIndices <= vValues) = vIndicesAbove(vGridIndices <= vValues) + 1;
vIndicesAbove(vIndicesAbove >= nGrid) = nGrid;
vGridIndicesAbove = vGrid(vIndicesAbove);

% Compute weights in the weighted sum
vWeightBelow = (vGridIndicesAbove - vValues) ./ (vGridIndicesAbove - vGridIndicesBelow);
vWeightBelow(vValues <= vGridIndicesBelow) = 1;
vWeightBelow(vValues >= vGridIndicesAbove) = 0;

vWeightAbove = (vValues - vGridIndicesBelow) ./ (vGridIndicesAbove - vGridIndicesBelow);
vWeightAbove(vValues <= vGridIndicesBelow) = 0;
vWeightAbove(vValues >= vGridIndicesAbove) = 1;
