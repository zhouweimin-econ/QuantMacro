function [hx,gx,F1,F2,F3,F4,par] = SGU_solver(F,mpar,par,p)
% This funtion solves for a competitive equilibrium defined as a zero of the
% function F which is written in Schmitt-Grohé Uribe form using the algorithm
% suggested  in Schmit-Grohé and Uribe (2004): "Solving dynamic general equilibrium models
% using a second-order approximation to the policy function"
%
%now do numerical linearization
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;
% tic
[Fb,~,~,~]  = F(State,State_m,Contr,Contr_m);%Fsys(State,State,Contr,Contr,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets);
% toc
%Use Schmitt Gohe Uribe Algorithm
% E[x' u']' =inv(A)*B*[x u]'
% A = [dF/dx' dF/du'], B =[dF/dx dF/du]
% A = [F1 F2]; B=[F3 F4]

F1 = zeros(mpar.numstates+mpar.numcontrols,mpar.numstates); %Tomorrows states do not affect error on controls and have unit effect on state error
F2 = zeros(mpar.numstates+mpar.numcontrols,mpar.numcontrols); %Jacobian wrt tomorrow's controls (TO BE FILLED)
F3 = zeros(mpar.numstates+mpar.numcontrols,mpar.numstates); % Jacobian wrt today's states (TO BE FILLED)
%Today's Value functions do not affect error on states and have unit effect
%on Value function error (LAST TWO COLUMNS TO BE FILLED: Aggregate Prices)
F4 = [zeros(mpar.numstates,mpar.numcontrols);eye(mpar.numcontrols,mpar.numcontrols)];


% disp('Use Schmitt Grohe Uribe Algorithm')
% disp(' A *E[xprime uprime] =B*[x u]')
% disp(' A = (dF/dxprimek dF/duprime), B =-(dF/dx dF/du)')
% A = [F1 F2]; B=[F3 F4]

% parameters which control numerical differentiation and Reiter solution
pnum     = p.NumWorkers;

packagesize = ceil(mpar.numstates/(3*pnum));
blocks      = ceil(mpar.numstates/packagesize);

% Absolute deviations
par.scaleval1 = 1e-5; %vector of numerical differentiation step sizes
par.scaleval2 = 1e-5; %vector of numerical differentiation step sizes

% jacobian wrt X', X
% tic
% disp('Computing Jacobian F1=DF/DXprime F3 =DF/DX')
% disp(['Total number of parallel blocks: ' num2str(blocks) '.'])
parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,mpar.numstates);
    DF1 = zeros(length(Fb),length(range));
    DF3 = zeros(length(Fb),length(range));
    cc = zeros(mpar.numcontrols,1);
    ss = zeros(mpar.numstates,1);
    for Xct = range
        X = zeros(mpar.numstates,1);
        h = par.scaleval1;
        X(Xct) = h;
        Fx = F(ss,X,cc,cc);%#ok Fsys(S_sp,X,C_sp,Cm_sp,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets);
        DF3(:,Xct-(bl-1)*packagesize)= (Fx-Fb)/h;
        Fx = F(X,ss,cc,cc);%Fsys(S_sp,Sm_sp,C_sp,Cm_sp,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets);
        DF1(:,Xct-(bl-1)*packagesize)= (Fx-Fb)/h;
    end    
    FF1{bl}=DF1;
    FF3{bl}=DF3;
%     disp(['Block number: ' num2str(bl) ' done.'])
end


for i=1:ceil(mpar.numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,mpar.numstates);
    F1(:,range)=FF1{i};
    F3(:,range)=FF3{i};
end
% toc

% jacobian wrt Y'
packagesize = ceil(mpar.numcontrols/(3*pnum));
blocks=ceil(mpar.numcontrols/packagesize);
% disp('Computing Jacobian F2 - DF/DYprime')
% disp(['Total number of parallel blocks: ' num2str(blocks) '.'])
parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,mpar.numcontrols);
    DF2=zeros(length(Fb),length(range));
    cc = zeros(mpar.numcontrols,1);
    ss = zeros(mpar.numstates,1);
    for Yct = range
        Y  = zeros(mpar.numcontrols,1);
        h = par.scaleval2;
        Y(Yct) = h;
        Fx = F(ss,ss,Y,cc);%#ok Fsys(NminusS_sp,X,C_sp,Cm_sp,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets);
        DF2(:,Yct-(bl-1)*packagesize)=(Fx-Fb)/h;
    end
    FF{bl}=DF2;
%     disp(['Block number: ' num2str(bl) ' done.'])
end

for i=1:ceil(mpar.numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,mpar.numcontrols);
    F2(:,range)=FF{i};
end
% toc
clear FF FF1 FF3;
% Derivative wrt Y (value functions today change only themselves)
%
cc = zeros(mpar.numcontrols,1);
ss = zeros(mpar.numstates,1);
for Yct = 0:mpar.oc-1
    Y  = zeros(mpar.numcontrols,1);
    h = par.scaleval2;
    Y(end-Yct) = h;
    Fx = F(ss,ss,cc,Y);%Fsys(S_sp,X,C_sp,Cm_sp,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets);
    F4(:,end-Yct)=(Fx-Fb)/h;
end

%% Schmidt-Grohé and Uribe code based on qz decomposition of the system
%  Code adapted from SGUs online ressources

[s, t, Q, Z] = qz(full([F1,F2]),full(-[F3, F4]));

relev = abs(diag(s))./abs(diag(t));
ll    = sort(relev);
slt   = relev>=1;
nk    = sum(slt); % Number of state Variables based on Eigenvalues
if nk>mpar.numstates
    if mpar.overrideEigen
        warning(['The Equilibrium is Locally Indeterminate, critical eigenvalue shifted to: ' num2str(ll(end-mpar.numstates))])
        slt = relev>ll(end-mpar.numstates);
        nk  = sum(slt);
    else
        error(['No Local Equilibrium Exists, last eigenvalue: ' num2str(ll(end-mpar.numstates))])
        %return
    end
elseif nk<mpar.numstates
    if mpar.overrideEigen
        warning(['No Local Equilibrium Exists, critical eigenvalue shifted to: ' num2str(ll(end-mpar.numstates))])
        slt = relev>ll(end-mpar.numstates);
        nk  = sum(slt);
    else
        error(['No Local Equilibrium Exists, last eigenvalue: ' num2str(ll(end-mpar.numstates))])
        %  return
    end
end

[s,t,~,Z] = ordqz(s,t,Q,Z,slt);

z21=Z(nk+1:end,1:nk);
z11=Z(1:nk,1:nk);
s11=s(1:nk,1:nk);
t11=t(1:nk,1:nk);

%Checks


if rank(z11)<nk
    warning('invertibility condition violated')
end
z11i=z11\eye(nk);
gx=real(z21*z11i);
hx=real(z11*(s11\t11)*z11i);

end
