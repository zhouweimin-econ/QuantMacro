function [P, y, x, MVndx, Moments0, p0] = tauchen(F, F0, Sigma, Ns, bandwidth)

% convert VAR(1) for y into Markov-Chain using Tauchen's method
% cool: Y can have arbitrary dimension!
% Y_t = F0 + F Y_{t-1} + e_t and E e_t e_t' = Sigma
% P is transition matrix
% y and x are midpoints of the Markov grid. x being the orthogonalized processes (using a choleski of Sigma)
% Ns is number of states *per variable* ! (Nstar = Ns^Ny)
% bandwidth is multiple of standard deviations which will be covered per variable
%
% Note: this function imposes same bandwidth and number of gridpoints per element of Y (makes algorithm more straightforward)

error(nargchk(3,5,nargin))

if nargin < 4
   Ns = 20; % number of states per variables
end

if nargin < 5
   bandwidth = 3; % number of standard deviations around average   
end

%% construct the range over which each of the compents may vary
% Y_t = F0 + F Y_{t-1} + Q \varepsilon_t
Ny       = size(F, 1);
Nstar    = Ns^Ny;
Q        = chol(Sigma)';
iQ       = inv(Q);
Iy       = eye(Ny);
EY       = inv(Iy - F) * F0;
%VarY     = disclyap(F, Q * Q');
VarY     = reshape(inv((eye(Ny^2) - kron(F,F)))*Sigma(:), Ny, Ny);

% X_t = Q^{-1} Y_t = F0x + Fx X_{t-1} + \varepsilon_t 
Fx       = iQ * F * Q;
F0x      = iQ * F0;
EX       = iQ * EY;
VarX     = iQ * VarY * iQ'; % = dlyap(F, Iy);
StdX     = sqrt(diag(VarX));

%% construct univariate grids for x (always midpoints!)
griddy  = repmat(NaN, Ns, Ny);
Ub      = EX + bandwidth * StdX;
Lb      = EX - bandwidth * StdX;
steps   = (Ub - Lb) / (Ns - 1);
for x = 1 : Ny
   griddy(:,x) = Lb(x) : steps(x) : Ub(x);
end

% griddy

%% Index for Multivariate Grid
MVndx                = gridMVndx(Ns, Ny);

% note: midpoints also used for conditional means! (see below)
x = griddy(MVndx);
y = x * Q';

endpoints = griddy(1:end-1,:) + repmat(steps' / 2, Ns-1,1);


%% conditional distributions
% note usage of griddy, not XuvgridMid!
condmean = x * Fx' + repmat(F0x',Nstar,1);
condstd  = ones(1, Ny);
P = repmat(NaN, Nstar, Nstar);
for s = 1 : Nstar
   probby   = diff([zeros(1,Ny); ...
      normcdf(endpoints, repmat(condmean(s,:), Ns-1, 1), repmat(condstd,Ns-1,1)); ...
      ones(1,Ny)]);
   P(s,:)   = prod(probby(MVndx), 2)';
end

%% construct unconditional distribution -- diagonalize VarX !!
if nargout > 5 
   if Ny > 1 
      warning('elmi:debu', 'p0 not correctly implemented! need MVgrid')
   end
   
   colly       = inv(chol(VarX)');
   uncondmean  = (colly * EX)';
   uncondstd   = sqrt(diag(colly * VarX * colly'))';
   ProbUV      = diff([zeros(1,Ny); ...
      normcdf(endpoints * colly', repmat(uncondmean, Ns-1, 1), repmat(uncondstd, Ns-1, 1)); ...
      ones(1,Ny)]);
   p0          = prod(ProbUV(MVndx), 2);
end

if nargout > 4
   Moments0.EY = EY;
   Moments0.EX = EX;
   Moments0.VarY = VarY;
   Moments0.VarX = VarX;
end


%% subfunction for Mulitvariate Index
function MVndx = gridMVndx(Ns, Ny)
% function MVndx = gridMVndx(Ns, Ny)
% to be used for constructing Xgrid  = griddy(MVndx);

Nstar = Ns^Ny;
MVndx = repmat(NaN, Nstar, Ny);
UVndx = (1:Ns)';
MVndx(:,1) = repmat(UVndx, Ns^(Ny-1), 1);

for j = 2 : Ny
   n1          = Ns ^ (j - 1);
   n2          = Ns ^ (Ny - j);
   MVndx(:,j)  = repmat(kron(UVndx, ones(n1,1)), n2, 1) + (j-1) * Ns;
end

