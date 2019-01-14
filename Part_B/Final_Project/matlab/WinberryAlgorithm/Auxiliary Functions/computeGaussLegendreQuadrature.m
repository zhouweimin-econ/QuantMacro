function [grid,weight] = computeGaussLegendreQuadrature(order)

% This function computes the notes and weights for the Gauss - Legendre quadrature
% of order "order." 
%
% Inputs
%   (1) order: order of the quadrature
%
% Outputs
%   (1) grid: nodes to evaluate the function on
%   (2) weights: weight in the approximation
%
% Jung Sakong, February 10, 2016

% Compute polynomial recursively
temp_poly = zeros(order+1,order+1);
temp_poly(1,1) = 1;
temp_poly(2,2) = 1;
for ii = 3:(order+1)
	temp_poly(ii,:) = (2*ii-3)/(ii-1)*[0,temp_poly(ii-1,1:end-1)]...
		- (ii-2)/(ii-1)*temp_poly(ii-2,:);
end
the_poly = fliplr(temp_poly(end,:)); % higher order coefficients first

% Solve for roots of the polynomial
grid = roots(the_poly);

% Compute weights 
temp_powers = zeros(order,order);
for ii = 1:order
	temp_powers(:,ii) = (order+1-ii)*grid.^(order-ii);
end
poly_prime = repmat(the_poly(1:end-1),[order 1]) ...
	.* temp_powers;
weight = 2 ./ ((1-grid.^2).*(sum(poly_prime,2)).^2);
