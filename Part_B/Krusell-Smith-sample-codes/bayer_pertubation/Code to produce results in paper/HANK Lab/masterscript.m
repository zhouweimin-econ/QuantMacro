% =========================================================================
% Part of the Matlab code to accompany the paper
% 'Solving heterogeneous agent models in discrete time with
% many idiosyncratic states by perturbation methods', CEPR Discussion Paper
% DP13071
% by Christian Bayer and Ralph Luetticke
% http://www.erc-hump-uni-bonn.net/
% =========================================================================
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The most recent version or successor of the above-mentioned paper is
% properly cited in all work and publications that benefit from insights
% derived from the Software.
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% =========================================================================

%% Initialize workspace and load directories
clear
clc
close all
restoredefaultpath

Computername='METHOD'
starttime=clock;
casenamemaster='SS_BASELINE_TWOASSET_TFP';
addpath(genpath('shared_functions'))
oldpath = path;

FindNewSS           = false;

%% Select aggregate shock
aggrshock           = 'TFP';
par.rhoS            = 0.75;     % Persistence
par.sigmaS          = 0.007;    % STD 

mpar.T(1,1)         = 1000;     % number of periods simulation
mpar.T(2,1)         = 0;        % burn-in phase
mpar.maxlag         = 40;       % number of periods IRF

%% Draw random seed for simulation across models
epsilon_master=randn(RandStream('mcg16807', 'Seed',20180621),mpar.T(1,1),1);

scontmaster       = 0;
for t=1:mpar.T(1,1)
    scontmaster(t+1) = scontmaster(t)*par.rhoS + par.sigmaS*epsilon_master(t);
end

%% Solve for Steady state
if FindNewSS
    disp('Solving Steady State by EGM')
    
    tic
    % Set parameters
    defineSS_pars
    
    mainskript_steadystate
    toc
else
    disp('Loading Steady State')
    
    load(casenamemaster)
    
end
%%
p = gcp('nocreate');
pool='local'; % select machine

if isempty(p)
    c = parcluster;
    parpool(pool,c.NumWorkers-1)
    p = gcp('nocreate');
end
%%
disp('Two Asset HANK - State Space Reduction')
cd('Two Asset HANK')
tic
mainskript
toc

%% Simulations
% DenHaan error
JDminus=joint_distr;
Brealized_RC=targets.B;
Krealized_RC=targets.K;
Hrealized_RC=par.H;
x=zeros(mpar.numstates,1);
RBminus=log(par.RB);
PIGamma= pinv(Gamma_state);

for t=1:mpar.T(1,1)-1
    ControlNOW=(gx*x)';
    xFOR=hx*x;
    ControlFOR=(gx*xFOR)';
    
    [JDminus,RBminus]  = denHaan_Reiter(JDminus,scontmaster(t),RBminus,...
        ControlFOR',  ControlNOW',...
        Yss,indexMUdct,indexVKdct,par,mpar,grid,targets,meshes,P_H,aggrshock);
    
    Brealized_RC(t+1)=meshes.m(:)'*JDminus(:);
    Krealized_RC(t+1)=meshes.k(:)'*JDminus(:);
    Hrealized_RC(t+1)=meshes.h(:)'*JDminus(:);
    % liquid assets
    aux_m = squeeze(sum(sum(JDminus,2),3));
    % illiquid assets
    aux_k =  squeeze(sum(sum(JDminus,1),3))' ;
    % human capital
    aux_h = squeeze(sum(sum(JDminus,1),2));
    Xaux=([aux_m(1:end);aux_k(1:end);aux_h(1:end)]-Xss(1:(mpar.nm+mpar.nk+mpar.nh)));
    
    x(1:end-2)=PIGamma*Xaux;
    x(end-1)=RBminus-log(par.RB);
    x(end)=scontmaster(t+1);
end
%% Produce IRFs and simulate
x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

plot_IRFs_prep


save(['../IRF_' casenamemaster],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

xlinear=x0;
IRF_state_sparse=[];
for t=1:mpar.T(1,1)
    xlinear(end)= scontmaster(t);
    IRF_state_sparse(:,t)=(MX*xlinear)';
    xlinear=hx*xlinear;
end

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-mpar.os,1:mpar.T(1,1));
Kaux=grid.k*IRF_distr(mpar.nm+(1:mpar.nk),:)+grid.K;
Baux=grid.m*IRF_distr(1:mpar.nm,:)+grid.B;

Kmat(:,1)=Kaux;
Imat(:,1)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));
DenHaanError(:,1)=100*log(Kaux./Krealized_RC);
DenHaanError(:,2)=100*log(Baux./Brealized_RC);

cd ..
% path(oldpath)

filename=['../latex/AvgAbsErrordenHaan_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f  \n', mean(abs(DenHaanError),1));
fclose(FID);
filename=['../latex/MaxAbsErrordenHaan_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f  \n',max(abs(DenHaanError),[],1));
fclose(FID);

plot_IRFs

%% 
disp('Solve Two Asset Model with Uncertainty Shocks')

casenamemaster='SS_BASELINE_TWOASSET_UNC';

%% Select aggregate shock
aggrshock           = 'Uncertainty';

par.rhoS            = 0.84;    % Persistence
par.sigmaS          = 0.54;    % STD 

%% SS
disp('Solving Steady State by EGM')

tic
% Set parameters
defineSS_pars_UNC

mainskript_steadystate
toc
%%
p = gcp('nocreate');
pool='local'; % select machine

if isempty(p)
    c = parcluster;
    parpool(pool,c.NumWorkers-1)
    p = gcp('nocreate');
end
%%
disp('Two Asset HANK - State Space Reduction')
cd('Two Asset HANK')
tic
mainskript
toc
%% Produce IRFs and simulate
x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

plot_IRFs_prep

save(['../IRF_' casenamemaster],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

cd ..

plot_IRFs
