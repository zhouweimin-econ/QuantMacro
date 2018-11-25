tic
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));

Xss=[squeeze(sum(sum(joint_distr,2),3)); ... % marginal distribution liquid
    squeeze(sum(sum(joint_distr,1),3))'; ... % marginal distribution illiquid
    squeeze(sum(sum(joint_distr,2),1)); ... % marginal distribution productivity
    log(par.RB); 0];

Yss=[invmutil(mutil_c(:)); invmutil(Vk(:)); log(par.Q); log(par.PI); log(Output);...
    log(par.G); log(par.W) ; log(par.R); log(par.PROFITS); log(par.N); log(targets.T);...
    log(grid.K); log(targets.B)];

%% Mapping for Histogram

Gamma_state=zeros(mpar.nm+mpar.nk+mpar.nh,mpar.nm+mpar.nk+mpar.nh-4);
for j=1:mpar.nm-1
    Gamma_state(1:mpar.nm,j)=-Xss(1:mpar.nm);
    Gamma_state(j,j)=1-Xss(j);
    Gamma_state(j,j)=Gamma_state(j,j) -sum(Gamma_state(1:mpar.nm,j));
end
bb=mpar.nm;
for j=1:mpar.nk-1
    Gamma_state(bb+(1:mpar.nk),bb+j-1)=-Xss(bb+(1:mpar.nk));
    Gamma_state(bb+j,bb-1+j)=1-Xss(bb+j);
    Gamma_state(bb+j,bb-1+j)=Gamma_state(bb+j,bb-1+j) -sum(Gamma_state(bb+(1:mpar.nk),bb-1+j));
end
bb=mpar.nm+mpar.nk;
for j=1:mpar.nh-2
    Gamma_state(bb+(1:mpar.nh-1),bb+j-2)=-Xss(bb+(1:mpar.nh-1));
    Gamma_state(bb+j,bb-2+j)=1-Xss(bb+j);
    Gamma_state(bb+j,bb-2+j)=Gamma_state(bb+j,bb-2+j) -sum(Gamma_state(bb+(1:mpar.nh-1),bb-2+j));
end

mpar.os = length(Xss) - (mpar.nm+mpar.nk+mpar.nh);
mpar.oc = length(Yss) - 2*(mpar.nm*mpar.nk*mpar.nh);

%%
[indexMUdct,DC1,DC2,DC3] = do_dct(invmutil(mutil_c(:)),mpar, 0.99999);

[indexVKdct] = do_dct(invmutil(Vk),mpar, 0.99999);

%%
aux= size(Gamma_state); %used for distributions

mpar.numstates   = aux(2)+mpar.os;
mpar.numcontrols = length(indexMUdct) + length(indexVKdct) + mpar.oc;
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;

