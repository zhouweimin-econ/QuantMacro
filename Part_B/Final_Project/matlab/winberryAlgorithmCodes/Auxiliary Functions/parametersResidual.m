function [value vDerivs] = parametersResidual(vParameters,mGridMoments);

% Computes the objective function and Jacobian for computing the parameters of the distribution, given
% moments to match
% 
% Inputs
%   (1) vParameters: candidate parameters of distribution (to be optimized over)
%	(2) mGridMoments: grid of centralized moments
%
% Outputs
%   (1) value: value of objective function
%	(2) vDerivs: Jacobian 
% 
% Thomas Winberry, October 12, 2015

global assetsMax assetsMin vQuadratureWeights nMeasure

% Value of function
value = vQuadratureWeights' * exp(mGridMoments * vParameters);

% Derivatives of function 
if nargout > 1

	vDerivs = sum(repmat(vQuadratureWeights,[1 nMeasure]) .* mGridMoments .* repmat(exp(mGridMoments * vParameters),...
		[1 nMeasure]),1)';

end

