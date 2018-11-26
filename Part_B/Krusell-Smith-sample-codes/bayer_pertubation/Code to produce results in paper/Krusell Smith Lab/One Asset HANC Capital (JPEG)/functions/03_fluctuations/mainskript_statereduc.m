
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));


Xss=[squeeze(sum(joint_distr,2)); ... % marginal distribution liquid
    squeeze(sum(joint_distr,1)'); ... % marginal distribution productivity
    0];

Yss=[invmutil(mutil_c(:)); log(par.Q); log(targets.Y);...
    log(par.W); log(par.R); log(par.N); log(grid.K)];

os = length(Xss) - (mpar.nm+mpar.nh);
oc = length(Yss) - (mpar.nm*mpar.nh);
%%
Gamma_state=zeros(mpar.nm+mpar.nh,mpar.nm+mpar.nh-2);
for j=1:mpar.nm-1
    Gamma_state(1:mpar.nm,j)=-Xss(1:mpar.nm);
    Gamma_state(j,j)=1-Xss(j);
    Gamma_state(j,j)=Gamma_state(j,j) -sum(Gamma_state(1:mpar.nm,j));
end
bb=mpar.nm;
for j=1:mpar.nh-1
    Gamma_state(bb+(1:mpar.nh),bb+j-1)=-Xss(bb+(1:mpar.nh));
    Gamma_state(bb+j,bb-1+j)=1-Xss(bb+j);
    Gamma_state(bb+j,bb-1+j)=Gamma_state(bb+j,bb-1+j) -sum(Gamma_state(bb+(1:mpar.nh),bb-1+j));
end


%%
aux = reshape(invmutil(mutil_c),[mpar.nm, mpar.nh]);
aux = mydct(aux,1); % do dct-transformation
aux = mydct(aux,2); % do dct-transformation
DC=aux(:);


[~,ind] = sort(abs(DC(:)),'descend');
DCaux=DC(:);
i = 1;
while norm(DCaux(ind(1:i)))/norm(DCaux) < 0.9999
    i = i + 1;
end
keepnumber = i;
indexMUdct=sort(ind(1:keepnumber));
%%
aux= size(Gamma_state); %used for distributions

mpar.numstates   = aux(2)+os;
mpar.numcontrols = length(indexMUdct)+oc;
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;
