function mPoly = computeChebyshev(nPower,vGrid);

% Computes Chebyshev polynomials up to order "nPower," evaluated along "vGrid"
% 
% Inputs
%   (1) nPower: order of polynomial
%	(2) vGrid: points along which to compute polynomials (must be column vector)
%
% Outputs
%   (1) mPoly: nGrid x nPower matrix of polynomials; entry (i,j) is the j_th order
%		Chebyshev polynomial evaluated at vGrid(i)
% 
% Thomas Winberry, January 19, 2016

% Compute grid size
[nGrid,~] = size(vGrid);

% Create polynomial
mPoly = ones(nGrid,nPower);
mPoly(:,2) = vGrid;
for iPower = 3:nPower
	mPoly(:,iPower) = 2 * vGrid .* mPoly(:,iPower-1) - mPoly(:,iPower-2);
end