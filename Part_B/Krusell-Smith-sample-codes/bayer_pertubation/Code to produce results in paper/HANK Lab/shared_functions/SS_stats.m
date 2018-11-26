% Steady state statistics
%%

[mesh.k,mesh.m]=meshgrid(grid.k,grid.m);
clear targets
targets.ShareBorrower=sum((grid.m<0).*sum(sum(joint_distr,2),3)');
targets.K=sum(grid.k.*sum(sum(joint_distr,1),3));
targets.B=grid.m*sum(sum(joint_distr,2),3);
grid.K=targets.K;
grid.B=targets.B;

JDredux = sum(joint_distr,3);
targets.BoverK     = targets.B/targets.K;

targets.L=grid.N*sum(grid.h'.*squeeze(sum(sum(joint_distr,1),2)));
targets.KY=targets.K/Output;
targets.BY=targets.B/Output;
targets.Y=Output;
BCaux_M=sum(sum(joint_distr,2),3);
targets.m_bc=BCaux_M(1,:);
targets.m_0=BCaux_M(grid.m==0);
BCaux_K=sum(sum(joint_distr,1),3);
targets.k_bc=BCaux_K(:,1);
aux_MK=sum(joint_distr,3);
targets.WtH_b0=sum(aux_MK(mesh.m==0 & mesh.k>0));
targets.WtH_bnonpos=sum(aux_MK(mesh.m<=0 & mesh.k>0));

targets.T =(1-par.tau).*W_fc.*grid.N +(1-par.tau).*Profits_fc;
par.G=targets.B*(1-par.RB/par.PI)+targets.T;
par.R=R_fc;
par.W=W_fc(1);
par.PROFITS=Profits_fc(1);
par.N=grid.N;
targets.GtoY=par.G/Output;

%% Ginis
% Net worth Gini
mplusk=mesh.k(:)*par.Q+mesh.m(:);
[mplusk, IX]       = sort(mplusk);
moneycapital_pdf   = JDredux(IX);
moneycapital_cdf   = cumsum(moneycapital_pdf);
targets.NegNetWorth= sum((mplusk<0).*moneycapital_pdf);

S                  = cumsum(moneycapital_pdf.*mplusk)';
S                  = [0 S];
targets.GiniW      = 1-(sum(moneycapital_pdf.*(S(1:end-1)+S(2:end))')/S(end));

% Liquid Gini
[liquid_sort, IX]  = sort(mesh.m(:));
liquid_pdf         = JDredux(IX);
liquid_cdf         = cumsum(liquid_pdf);
targets.Negliquid  = sum((liquid_sort<0).*liquid_pdf);

S                  = cumsum(liquid_pdf.*liquid_sort)';
S                  = [0 S];
targets.GiniLI      = 1-(sum(liquid_pdf.*(S(1:end-1)+S(2:end))')/S(end));

% Illiquid Gini
[illiquid_sort, IX] = sort(mesh.k(:));
illiquid_pdf        = JDredux(IX);
illiquid_cdf        = cumsum(illiquid_pdf);
targets.Negliquid   = sum((illiquid_sort<0).*illiquid_pdf);

S                   = cumsum(illiquid_pdf.*illiquid_sort)';
S                   = [0 S];
targets.GiniIL      = 1-(sum(illiquid_pdf.*(S(1:end-1)+S(2:end))')/S(end));

%%   MPCs
[meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

NW=par.gamma/(1+par.gamma).*(par.N/par.H).*par.W;
WW=NW*ones(mpar.nm,mpar.nk,mpar.nh); %Wages
WW(:,:,end)=par.PROFITS*par.profitshare;
% MPC
WW_h=squeeze(WW(1,1,:));
WW_h_mesh=squeeze(WW(:,:,:).*meshes.h);

grid_h_aux=grid.h;

MPC_a_m = zeros(mpar.nm,mpar.nk,mpar.nh);
MPC_n_m = zeros(mpar.nm,mpar.nk,mpar.nh);
for kk=1:mpar.nk
    for hh=1:mpar.nh
        MPC_a_m(:,kk,hh)=gradient(squeeze(c_a_guess(:,kk,hh)))./gradient(grid.m)';
        MPC_n_m(:,kk,hh)=gradient(squeeze(c_n_guess(:,kk,hh)))./gradient(grid.m)';
    end
end

MPC_a_m=MPC_a_m.*(WW_h_mesh./c_a_guess);
MPC_n_m=MPC_n_m.*(WW_h_mesh./c_n_guess);

MPC_a_h = zeros(mpar.nm,mpar.nk,mpar.nh);
MPC_n_h = zeros(mpar.nm,mpar.nk,mpar.nh);
for mm=1:mpar.nm
    for kk=1:mpar.nk
        MPC_a_h(mm,kk,:)=gradient(squeeze(log(c_a_guess(mm,kk,:))))./gradient(log(WW_h'.*grid_h_aux))';
        MPC_n_h(mm,kk,:)=gradient(squeeze(log(c_n_guess(mm,kk,:))))./gradient(log(WW_h'.*grid_h_aux))';
    end
end

EMPC_h=joint_distr(:)'*(par.nu.*MPC_a_h(:)+(1-par.nu).*MPC_n_h(:));
EMPC_m=joint_distr(:)'*(par.nu.*MPC_a_m(:)+(1-par.nu).*MPC_n_m(:));

EMPC_a_h=joint_distr(:)'*MPC_a_h(:);
EMPC_a_m=joint_distr(:)'*MPC_a_m(:);

EMPC_n_h=joint_distr(:)'*MPC_n_h(:);
EMPC_n_m=joint_distr(:)'*MPC_n_m(:);

targets.Insurance_coeff=[1-EMPC_h 1-EMPC_m;
    1-EMPC_a_h 1-EMPC_a_m;
    1-EMPC_n_h 1-EMPC_n_m];
