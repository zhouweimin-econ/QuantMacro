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
casenamemaster='SS_BASELINE_HANC';
addpath(genpath('shared_functions'))
oldpath = path;

p = gcp('nocreate');
pool='local'; % select machine

if isempty(p)
    c = parcluster;
    parpool(pool,c.NumWorkers)
    p = gcp('nocreate');
end

%% Select aggregate shock
aggrshock           = 'TFP';
par.rhoS            = 0.75;     % Persistence 
par.sigmaS          = 0.007;    % STD 
mpar.ns             = 2;        % number of TFP states
mpar.T(1,1)         = 1500;     % number of periods
mpar.T(2,1)         = 500;      % burn-in phase
mpar.maxlag         = 40;       % number of periods IRF

mpar.nK             = 3;        % number of grid points for KS


%% Draw random seed for simulation across models
pr_s_master=rand(RandStream('mcg16807', 'Seed',20180621),1,mpar.T(1,1));

[P_S, grid.s] = rouwen(par.rhoS, 0, par.sigmaS/sqrt(1-par.rhoS^2), mpar.ns);
grid.s=exp(grid.s');
PS_S=cumsum(P_S,2);

smaster(1)        =ceil(mpar.ns/2);
for t=1:mpar.T(1,1)   
    smaster(t+1) = min(mpar.ns,sum(PS_S(smaster(t),:)<pr_s_master(t))+1);
end

%% Solve for Steady state
disp('Solving Steady State by EGM')
tic
% Set parameters
defineSS_pars

mainskript_steadystate
toc
%%
disp('Reiter method with state space reduction')
cd('One Asset HANC Capital (JPEG)')
tic
mainskript
toc

%% Simulations
% DenHaan error
JDminus=joint_distr;
Krealized_RC=targets.K;
Hrealized_RC=par.H;
x=zeros(mpar.numstates,1);

PIGamma= pinv(Gamma_state);

for t=1:mpar.T(1,1)-1
    ControlNOW=(gx*x)';
    xFOR=hx*x;
    ControlFOR=(gx*xFOR)';
    
    [JDminus]  = denHaan_Reiter(JDminus,log(grid.s(smaster(t))),...
        ControlFOR',  ControlNOW',...
        Yss,indexMUdct,par,mpar,grid,meshes,P_H,aggrshock,oc);

    Krealized_RC(t+1)=meshes.m(:)'*JDminus(:);
    Hrealized_RC(t+1)=meshes.h(:)'*JDminus(:);
    % liquid assets
    aux_m = squeeze(sum(JDminus,2));
    % human capital
    aux_h = squeeze(sum(JDminus,1))';
    Xaux=([aux_m(1:end);aux_h(1:end)]-Xss(1:(mpar.nm+mpar.nh)));
    
    x(1:end-1)=PIGamma*Xaux;
    x(end)=log(grid.s(smaster(t+1)));    
end
cd ..
path(oldpath)
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

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-1,1:mpar.maxlag);

IRF_H=100*grid.h(1:end)*IRF_distr(mpar.nm+(1:mpar.nh),2:end)/par.H;
K=grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I=100*(I/(par.delta*grid.K)-1);
IRF_K=100*(K(1:end)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);

% plot_IRFs

save(['IRF_' casename],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

xlinear=x0;
IRF_state_sparse=[];
for t=1:mpar.T(1,1)
    xlinear(end)= log(grid.s(smaster(t)));
    IRF_state_sparse(:,t)=(MX*xlinear)';
    xlinear=hx*xlinear;
end

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-1,1:mpar.T(1,1));

Kaux=(grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K)';
Kmat(:,1)=Kaux;
Imat(:,1)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));

DenHaanError(:,1)=100*log(Kaux./Krealized_RC');

%%
disp('Reiter method without reduction')
cd('One Asset HANC Capital (FULL)')
tic
mainskript
toc

%% Simulations
% DenHaan error
JDminus=joint_distr;
Krealized_RF=targets.K;
Hrealized_RF=par.H;
x=zeros(mpar.numstates,1);

PIGamma= pinv(Gamma_state);

for t=1:mpar.T(1,1)-1
    ControlNOW=(gx*x)';
    xFOR=hx*x;
    ControlFOR=(gx*xFOR)';
    
    [JDminus]  = denHaan_Reiter(JDminus,log(grid.s(smaster(t))),...
        ControlFOR',  ControlNOW',...
        Yss,indexMUdct,par,mpar,grid,meshes,P_H,aggrshock,oc);

    Krealized_RF(t+1)=meshes.m(:)'*JDminus(:);
    Hrealized_RF(t+1)=meshes.h(:)'*JDminus(:);

    Xaux=([JDminus(:)]-Xss(1:(mpar.nm*mpar.nh)));
    
    x(1:end-1)=PIGamma*Xaux;
    x(end)=log(grid.s(smaster(t+1)));    
end
cd ..
path(oldpath)
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


IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-os,1:mpar.maxlag);

IRF_distr=reshape(IRF_distr,[mpar.nm mpar.nh mpar.maxlag]);
IRF_distr_m=squeeze(sum(IRF_distr,2));
IRF_distr_h=squeeze(sum(IRF_distr,1));
IRF_H=100*grid.h(1:end)*IRF_distr_h(1:end,2:end)/par.H;
K=grid.m*IRF_distr_m(:,1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I=100*(I/(par.delta*grid.K)-1);
IRF_K=100*(K(1:end)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);

% plot_IRFs

save(['IRF_' casename],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

xlinear=x0;
IRF_state_sparse=[];
for t=1:mpar.T(1,1)
    xlinear(end)= log(grid.s(smaster(t)));
    IRF_state_sparse(:,t)=(MX*xlinear)';
    xlinear=hx*xlinear;
end

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-os,1:mpar.T(1,1));

IRF_distr=reshape(IRF_distr,[mpar.nm mpar.nh mpar.T(1,1)]);
IRF_distr_m=squeeze(sum(IRF_distr,2));
IRF_distr_h=squeeze(sum(IRF_distr,1));
Kaux=grid.m*IRF_distr_m(:,1:end)+grid.K;
Kmat(:,2)=Kaux;
Imat(:,2)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));

DenHaanError(:,2)=100*log(Kaux./Krealized_RF);

%%
disp('Krusell Smith method')
cd('One Asset HANC Capital (KS)')
tic
mainskript
toc
cd ..
path(oldpath)
addpath(genpath('One Asset HANC Capital (KS)/functions'))

[IRF_K_KS, IRF_S_KS]=KSsimulationIRF(k_star,P_H,joint_distr,grid,mesh,mpar,par);

[Kmat(:,3)]=KSsimulation(k_star,P_H,PS_S,pr_s_master(:),joint_distr,grid,mesh,mpar);
Kaux=Kmat(:,3);
Imat(:,3)=(Kaux(2:end)-(1-par.delta)*Kaux(1:end-1));

Klinear_KS=zeros(mpar.T(1,1),1);
Klinear_KS(1)=targets.K;
for t=1:mpar.T(1,1)-1
    
    Klinear_KS(t+1)=exp(par.beta_K(1,smaster(t))   + par.beta_K(2,smaster(t)).*log(Klinear_KS(t)));
%             Klinear_KS(t+1)     = exp(par.beta_K(1)   + par.beta_K(2)*log(grid.s(smaster(t))) + par.beta_K(3)*log(Klinear_KS(t)));

    
end

DenHaanError(:,3)=100*log(Kaux./Klinear_KS);


path(oldpath)
%% Save simulation results
figurename=['K_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:mpar.T(1,1)-mpar.T(2,1),Kmat(mpar.T(2,1)+1:end,1),'--','LineWidth',3.5)
hold on
plot(1:mpar.T(1,1)-mpar.T(2,1),Kmat(mpar.T(2,1)+1:end,2),'-o','LineWidth',3.5)
plot(1:mpar.T(1,1)-mpar.T(2,1),Kmat(mpar.T(2,1)+1:end,3),'-','LineWidth',3.5)
legend('Reiter-Reduction','Reiter-Full','Krusell-Smith')
ylabel('Capital','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.1 0.15 0.85 0.825])
printpdf(gcf,['../latex/' figurename])

figurename=['I_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:mpar.T(1,1)-1,Imat,'LineWidth',3.5)
legend('Reiter-Reduction','Reiter-Full','Krusell-Smith')
ylabel('Investment','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.15 0.15 0.8 0.8])
printpdf(gcf,['../latex/' figurename])

Kmat_normalized=Kmat./repmat(mean(Kmat(mpar.T(2,1):end,:),1),[mpar.T(1,1) 1]);
RCtoKS=100*log(Kmat_normalized(mpar.T(2,1):end,1)./Kmat_normalized(mpar.T(2,1):end,3));
RFtoKS=100*log(Kmat_normalized(mpar.T(2,1):end,2)./Kmat_normalized(mpar.T(2,1):end,3));

ErrorTabletoKSNorm(1,:)=mean(abs([RCtoKS RFtoKS]));
ErrorTabletoKSNorm(2,:)=max(abs([RCtoKS RFtoKS]));

filename=['../latex/AvgAbsErrortoKSNorm_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f \n', ErrorTabletoKSNorm(1,1:2));
fclose(FID);
filename=['../latex/MaxAbsErrortoKSNorm_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f \n', ErrorTabletoKSNorm(2,1:2));
fclose(FID);

figurename=['K_Normalized_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:mpar.T(1,1),Kmat_normalized(:,1),'--','LineWidth',3.5)
hold on
plot(1:mpar.T(1,1),Kmat_normalized(:,2),'-o','LineWidth',3.5)
plot(1:mpar.T(1,1),Kmat_normalized(:,3),'-','LineWidth',3.5)
legend('Reiter-Reduction','Reiter-Full','Krusell-Smith')
ylabel('Percent deviations','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.15 0.15 0.8 0.8])
printpdf(gcf,['../latex/' figurename])

figurename=['K_Normalized_Comparison_Short_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:100,Kmat_normalized(mpar.T(2,1)+1:mpar.T(2,1)+100,1),'--','LineWidth',3.5)
hold on
plot(1:100,Kmat_normalized(mpar.T(2,1)+1:mpar.T(2,1)+100,2),'o','LineWidth',3.5)
plot(1:100,Kmat_normalized(mpar.T(2,1)+1:mpar.T(2,1)+100,3),'-','LineWidth',3.5)
legend('Reiter-Reduction','Reiter-Full','Krusell-Smith')
ylabel('Percent deviations','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.15 0.15 0.8 0.8])
printpdf(gcf,['../latex/' figurename])

RCtoKS=100*log(Kmat(:,1)./Kmat(:,3));
RFtoKS=100*log(Kmat(:,2)./Kmat(:,3));

ErrorTabletoKS(1,:)=mean(abs([RCtoKS RFtoKS]));
ErrorTabletoKS(2,:)=max(abs([RCtoKS RFtoKS]));


filename=['../latex/AvgAbsErrortoKS_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f \n', ErrorTabletoKS(1,1:2));
fclose(FID);
filename=['../latex/MaxAbsErrortoKS_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f \n', ErrorTabletoKS(2,1:2));
fclose(FID);

RCtoRF=100*log(Kmat(:,1)./Kmat(:,2));
ErrorTabletoRF(1,:)=mean(abs([RCtoRF]));
ErrorTabletoRF(2,:)=max(abs([RCtoRF]));

filename=['../latex/AvgAbsErrortoRF_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f \n', ErrorTabletoRF(1,:));
fclose(FID);
filename=['../latex/MaxAbsErrortoRF_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f \n', ErrorTabletoRF(2,:));
fclose(FID);

filename=['../latex/AvgAbsErrordenHaan_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f & %8.4f \n', mean(abs(DenHaanError),1));
fclose(FID);
filename=['../latex/MaxAbsErrordenHaan_' casenamemaster '_' aggrshock '.tex'];
FID = fopen(filename, 'w');
fprintf(FID, '   %8.4f & %8.4f & %8.4f \n',max(abs(DenHaanError),[],1));
fclose(FID);
%% Save IRFs

IRF_K_KSaux=100*(IRF_K_KS(mpar.T(2,1)+1:end)./IRF_K_KS(mpar.T(2,1))-1);
Iaux=(IRF_K_KS(2:end)-(1-par.delta)*IRF_K_KS(1:end-1));
IRF_I_KSaux=100*(Iaux(mpar.T(2,1)+1:end)./Iaux(mpar.T(2,1))-1);

load IRF_SS_BASELINE_HANC_JPEG
IRF_K_RC=IRF_K';
IRF_I_RC=IRF_I';

load IRF_SS_BASELINE_HANC_FULL
IRF_K_RF=IRF_K';
IRF_I_RF=IRF_I';

IRF_K=[IRF_K_RC IRF_K_RF IRF_K_KSaux];

IRF_I=[IRF_I_RC IRF_I_RF IRF_I_KSaux];


figurename=['IRF_K_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:mpar.maxlag,IRF_K(:,1),'--','LineWidth',3.5)
hold on
plot(1:mpar.maxlag,IRF_K(:,2),'o','LineWidth',3.5)
plot(1:mpar.maxlag,IRF_K(:,3),'-','LineWidth',3.5)
ylim([0 0.25])
legend('Reiter-Reduction','Reiter-Full','Krusell-Smith','Location','SouthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.125 0.15 0.85 0.825])
printpdf(gcf,['../latex/' figurename])

figurename=['IRF_I_Comparison_' casenamemaster  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 1500 1000])
plot(1:mpar.maxlag-1,IRF_I(:,1),'--','LineWidth',3.5)
hold on
plot(1:mpar.maxlag-1,IRF_I(:,2),'o','LineWidth',3.5)
plot(1:mpar.maxlag-1,IRF_I(:,3),'-','LineWidth',3.5)

legend('Reiter-Reduction','Reiter-Full','Krusell-Smith','Location','NorthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.125 0.15 0.85 0.825])
printpdf(gcf,['../latex/' figurename])

%% Produce results for uncertainty shock
disp('Solve KS model with uncertainty shocks')

casenamemaster='SS_BASELINE_HANC_UNC';

% Select aggregate shock
aggrshock           = 'Uncertainty';
par.rhoS            = 0.84;    % Persistence of variance
par.sigmaS          = 0.54;    % STD of variance shocks

%% Solve for Steady state
disp('Solving Steady State by EGM')
tic
% Set parameters
defineSS_pars_UNC

mainskript_steadystate
toc
%%
disp('Reiter method with State Space Reduction')
cd('One Asset HANC Capital (JPEG)')
tic
mainskript
toc
cd ..
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

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-1,1:mpar.maxlag);

IRF_H=100*grid.h(1:end)*IRF_distr(mpar.nm+(1:mpar.nh),2:end)/par.H;
K=grid.m*IRF_distr((1:mpar.nm),1:end)+grid.K;
I=(K(2:end)-(1-par.delta)*K(1:end-1));
IRF_I=100*(I/(par.delta*grid.K)-1);
IRF_K=100*(K(1:end)./grid.K-1);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+1,1:end-1);
Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_C=100*((Y-I)./(targets.Y-par.delta*grid.K)-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_Q=100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_R=100*IRF_state_sparse(end-oc+4,1:end-1);

casename=casenamemaster;

save(['IRF_' casename],'IRF_K','IRF_I','IRF_S','IRF_Y','IRF_C','IRF_Q','IRF_N','IRF_R','IRF_H')

plot_IRFs
