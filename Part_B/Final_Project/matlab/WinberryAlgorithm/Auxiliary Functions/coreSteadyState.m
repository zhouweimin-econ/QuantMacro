% Computes market clearing capital stock and associated distribution and 
% decision rules in steady state
%
% Thomas Winberry, July 26th, 2016

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
computeGrids;

% Polynomials over grids (only if using polynomials to approximate conditional expectation)
if splineOpt == 0
	computePolynomials;
end

%----------------------------------------------------------------
% Compute initial guess of market-clearing capital stock using
% histogram approximation of distribution, from Young (2010)
%----------------------------------------------------------------

t0 = tic;
fprintf('Computing initial guess from histogram...\n')

% Solve for market clearing capital stock
f = @(capital) computeMCResidualHistogram(capital);
options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-2); % In older versions of MATLAB, use: options = optimset('Display',displayOpt); 
[aggregateCapitalInit,err,exitflag] = fsolve(f,1.01*kRepSS,options);

% Return exitflag if market clearing not solved
if exitflag < 1
    check = 1;
    return; 
end	

aggregateCapital = aggregateCapitalInit;
if strcmp(displayOpt,'iter-detailed') == 1
    fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
end

%----------------------------------------------------------------
% Compute moments of histogram to use as initial guess for parametric family
%----------------------------------------------------------------

% Compute histogram
[~, mHistogram] = computeMCResidualHistogram(aggregateCapital);

% Compute moments from histogram
mMomentsHistogram = zeros(nEpsilon,nMeasure);
aGridMoments = zeros(nEpsilon,nAssetsQuadrature,nMeasure); % grid for computing PDF

for iEpsilon = 1 : nEpsilon
	
	% First moment (uncentered)
	mMomentsHistogram(iEpsilon,1) = sum(vAssetsGridFine' .* (mHistogram(iEpsilon,:) ./ ...
		sum(mHistogram(iEpsilon,:))));
	aGridMoments(iEpsilon,:,1) = vAssetsGridQuadrature - mMomentsHistogram(iEpsilon,1);
		
	% Higher order moments (centered)
	for iMoment = 2 : nMeasure
		mMomentsHistogram(iEpsilon,iMoment) = sum(((vAssetsGridFine' - mMomentsHistogram(iEpsilon,1)) .^ iMoment) .* ...
			(mHistogram(iEpsilon,:) ./ sum(mHistogram(iEpsilon,:))));
		aGridMoments(iEpsilon,:,iMoment) = (vAssetsGridQuadrature' - mMomentsHistogram(iEpsilon,1)) .^ ...
			iMoment - mMomentsHistogram(iEpsilon,iMoment);
	end	
	
end

% Mass at borrowing constraint
mHatHistogram = [mHistogram(1,1) / sum(mHistogram(1,:));mHistogram(2,1) / sum(mHistogram(2,:))];

%----------------------------------------------------------------
% Compute market-clearing capital stock from parametric family
%----------------------------------------------------------------

t0 = tic; 
fprintf('Compute steady state from parametric family...\n')

% Solve for market clearing capital stock
f = @(capital) computeMCResidualPolynomials(capital,mMomentsHistogram,aGridMoments,mHatHistogram);
options = optimoptions('fsolve','Display',displayOpt,'TolFun',1e-2);
if abs(f(aggregateCapitalInit)) > 1e-4
	[aggregateCapital,err,exitflag] = fsolve(f,aggregateCapitalInit,options);
end

% Return error if market clearing not solved
if exitflag < 1
    check = 1;
    return; 
end	

if strcmp(displayOpt,'iter-detailed') == 1
    fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
end

%----------------------------------------------------------------
% Compute other objects from steady state
%----------------------------------------------------------------

[~,mCoefficients,mParameters,mMoments,mHat] = ...
    computeMCResidualPolynomials(aggregateCapital,mMomentsHistogram,aGridMoments,mHatHistogram);