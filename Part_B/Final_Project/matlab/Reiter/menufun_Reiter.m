%{
 ________________________________________________________
    
 Final Project of Part B Quantitative Macroeconomics, UAB
 
 Weimin Zhou
 Jan, 2019
 ________________________________________________________
 
 Descriptions:
  the main file
    Heterogenous Agent DSGE model using Reiter Method
    Model set-up: Krusell-Smith Model
    codes references see my report file: Final_Project.pdf
  Algorithm:
  Step.1 compute the steady states with a sequence of guess and parameters
  Step.2 solve for the model dynamics (using 1st order perturbation)
  Step.3 calculate impulse responses
%}

clear all; clc;
global n indx nodes Grid
global betta alf sig delta Zbar rhoz sigma_epsz TransMat_IdShock amin amax
global coeffguess distrguess

% set up KS parameters
parameters;
%% Step.1  compute the steady state

% initial guess
K_stst = 1;
C_stst = (Zbar*K_stst^alf-delta*K_stst);
Inv_stst = delta*K_stst;
Y_stst = Zbar*K_stst^alf;
W_stst = (1-alf)*Y_stst;
R_stst = 1+alf*Y_stst/K_stst-delta;

bettaguess = log(0.98);
coeffguess = log((R_stst-1)*Grid(:,1)+W_stst*Grid(:,2));
distrguess = 1/n.Agents * ones(n.Agents,1);
test=model_stst(bettaguess);

% solve for steady steate
opt_options = optimset('Display','iter');
betta_sol = fzero(@model_stst,bettaguess,opt_options);
[~,coeff_stst,distr_stst] = model_stst(betta_sol); % recover other elements of solution

betta = exp(betta_sol);

stst = [coeff_stst; distr_stst; ...
        log([C_stst; Inv_stst; K_stst; R_stst; Y_stst; 1])];

%% Step.2 solve for the model dynamics (using 1st order perturbation)
% i.e. obtain the state-space: x(t) = G1*x(t-1) + impact*eps(t);
doAutDiff = 0;  % choose whether to use automatic differentiation
[G1,C,impact] = linear_solution(@model_equations,stst,indx,doAutDiff);

%% Step.3 calculate impulse responses
horizonmax   = 40;
scale_shocks = ones(n.Shocks,1);
scale_y      = ones(n.Var,1);
Resp_mat     = impresp_sims(G1,impact,horizonmax,scale_shocks,scale_y);

% store response of distribution
Resp_distr = distr_stst + Resp_mat(indx.current(n.Agents+1:2*n.Agents),:)/100;

% retain only rows of aggregate variables
Resp_mat = Resp_mat(indx.AggVar(1:n.Aggr),:);

%% Report the results

% Plot steady state distribution
figure(1)
area(nodes.assets, sum(reshape(distr_stst,n.I,n.S),2));
xlabel('Assets'); 
title('Steady State Distribution');

% Plot impulse responses of aggregate 
time=0:size(Resp_mat,2)-1;
nVarIRF = size(Resp_mat,1);
set(0,'DefaultAxesGridLineStyle',':');

for ishock=1:n.Shocks
    figure(100+ishock)
    for ivar=1:nVarIRF-n.Shocks
        subplot(ceil((nVarIRF-n.Shocks+1)/2),2,ivar);
        plot(time(:),Resp_mat(ivar,:,ishock),'o-','MarkerSize',5,'LineWidth',1.5);
        xlabel('time');
        title(vlab(ivar,:));
        axis([-0.5 horizonmax -Inf Inf]);
        axis 'auto y';
    end
    subplot(ceil((nVarIRF-n.Shocks+1)/2),2,ivar+1);
    plot(time(:),Resp_mat(ivar+ishock,:,ishock),'o-','MarkerSize',5,'LineWidth',1.5);
    xlabel('time');
    title(vlab(ivar+ishock,:));
    axis([-0.5 horizonmax -Inf Inf]);
    axis 'auto y';
end