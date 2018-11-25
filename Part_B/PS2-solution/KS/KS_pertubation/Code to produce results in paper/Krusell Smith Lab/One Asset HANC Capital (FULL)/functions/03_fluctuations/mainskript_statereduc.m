
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));


Xss=[joint_distr(:); 0];

Yss=[invmutil(mutil_c(:)); log(par.Q); log(targets.Y);...
     log(par.W); log(par.R); log(par.N); log(grid.K)];
 
os = length(Xss) - (mpar.nm*mpar.nh);
oc = length(Yss) - (mpar.nm*mpar.nh);
%%
Gamma_state=zeros(mpar.nm*mpar.nh,mpar.nm*mpar.nh-1);
for j=1:mpar.nm*mpar.nh-1
    Gamma_state(1:mpar.nm*mpar.nh,j)=-Xss(1:mpar.nm*mpar.nh);
    Gamma_state(j,j)=1-Xss(j);
    Gamma_state(j,j)=Gamma_state(j,j) -sum(Gamma_state(1:mpar.nm*mpar.nh,j));
end

%%
    indexMUdct=1:mpar.nm*mpar.nh;
%%
aux= size(Gamma_state); %used for distributions

mpar.numstates   = aux(2)+os;
mpar.numcontrols = length(indexMUdct)+oc;
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;
