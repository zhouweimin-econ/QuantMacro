% Computes and analyzes steady state with no aggregate shocks
%
% Thomas Winberry, July 26th, 2016

clear all
close all
clc

oldFolder = cd('./Auxiliary Functions');

%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Compute Steady State
%----------------------------------------------------------------

% Solve for steady state capital stock, distribution, and decision rules
coreSteadyState;

% Compute decision rules along fine grid for analysis
[~,mHistogram,mAssetsPrime,mConsumption] = computeMCResidualHistogram(aggregateCapital);

% Compute density along fine grid
mDistributionFine = zeros(nEpsilon,nAssetsFine);
for iEpsilon = 1 : nEpsilon

	% First moment (uncentered)
	mGridMoments = zeros(nAssetsFine,nMeasure);
	mGridMoments(:,1) = (vAssetsGridFine - mMoments(iEpsilon,1));
	
	% Higher order moments (centered)
	for iMoment = 2 : nMeasure
		mGridMoments(:,iMoment) = (vAssetsGridFine - mMoments(iEpsilon,1)) .^ iMoment - ...
			mMoments(iEpsilon,iMoment);
	end
	
	% Compute density away from borrowing constraint
	mDistributionFine(iEpsilon,:) = mParameters(iEpsilon,1) * exp(mGridMoments * ...
		mParameters(iEpsilon,2:nMeasure+1)');
		
	% Mass at borrowing constraint
	%mDistributionFine(iEpsilon,1) = mHat(iEpsilon,1); % Commented out for now; need fine quadrature grid to capture correctly
		
end

%----------------------------------------------------------------
% Plot results 
%----------------------------------------------------------------

% Savings function
figure
hold on
plot(vAssetsGridFine,mAssetsPrime(1,:),'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vAssetsGridFine,mAssetsPrime(2,:),'linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vAssetsGridFine,vAssetsGridFine,'k--','linewidth',1)
xlabel('Assets, $a$','interpreter','latex')
ylabel('Savings, $s(\varepsilon,a)$','interpreter','latex')
xlim([aaBar .9*assetsMax])
title('Savings Decision Rule')
legend('Unemployed','Employed','location','southeast')
grid on
set(gcf,'color','w')
hold off

% Consumption function
figure
hold on
plot(vAssetsGridFine,mConsumption(1,:),'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vAssetsGridFine,mConsumption(2,:),'linewidth',1.5,'color',[8/255,62/255,118/255])
xlabel('Assets, $a$','interpreter','latex')
ylabel('Consumption, $c(\varepsilon,a)$','interpreter','latex')
xlim([aaBar .9*assetsMax])
title('Consumption Decision rule')
legend('Unemployed','Employed','location','southeast')
grid on
set(gcf,'color','w')
hold off

% Distribution of unemployed households
figure
hold on
plot(vAssetsGridFine,mHistogram(1,:) / sum(mHistogram(1,:)),...
	'linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vAssetsGridFine,mDistributionFine(1,:) ./ sum(mDistributionFine(1,:)),...
	'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
xlabel('Assets, $a$','interpreter','latex')
ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
xlim([aaBar .9*assetsMax])
title('Invariant Distribution of Households (Unemployed)')
legend('Histogram','Parametric Family','location','northeast')
grid on
set(gcf,'color','w')
hold off

% Distribution of employed households
figure
hold on
plot(vAssetsGridFine,mHistogram(2,:) / sum(mHistogram(2,:)),...
	'linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vAssetsGridFine,mDistributionFine(2,:) ./ sum(mDistributionFine(2,:)),...
	'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
xlabel('Assets, $a$','interpreter','latex')
ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex')
xlim([aaBar .9*assetsMax])
title('Invariant Distribution of Households (Employed)')
legend('Histogram','Parametric Family','location','northeast')
grid on
set(gcf,'color','w')
hold off

cd(oldFolder)