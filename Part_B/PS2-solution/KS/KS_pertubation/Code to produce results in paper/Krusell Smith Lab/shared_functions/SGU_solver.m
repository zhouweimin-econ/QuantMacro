function [hx,gx,F1,F2,F3,F4,par] = SGU_solver(F, numstates,numcontrols,oc,mpar,par,p)
% This funtion solves for a competitive equilibrium defined as a zero of the
% function F which is written in Schmitt-Grohé Uribe form using the algorithm
% suggested  in Schmit-Grohé and Uribe (2004): "Solving dynamic general equilibrium models
% using a second-order approximation to the policy function"
%
%now do numerical linearization
State       = zeros(numstates,1);
State_m     = State;
Contr       = zeros(numcontrols,1);
Contr_m     = Contr;

[Fb,~,~,~]  = F(State,State_m,Contr,Contr_m);%Fsys(State,State,Contr,Contr,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets);

%Use Schmitt Gohe Uribe Algorithm
% E[x' u']' =inv(A)*B*[x u]'
% A = [dF/dx' dF/du'], B =[dF/dx dF/du]
% A = [F1 F2]; B=[F3 F4]

% p = gcp('nocreate');
F1 = zeros(numstates+numcontrols,numstates); %Tomorrows states do not affect error on controls and have unit effect on state error
F2 = zeros(numstates+numcontrols,numcontrols); %Jacobian wrt tomorrow's controls (TO BE FILLED)
F3 = zeros(numstates+numcontrols,numstates); % Jacobian wrt today's states (TO BE FILLED)
%Today's Value functions do not affect error on states and have unit effect
%on Value function error (LAST TWO COLUMNS TO BE FILLED: Aggregate Prices)
F4 = [zeros(numstates,numcontrols);eye(numcontrols,numcontrols)];


% disp('Use Schmitt Grohe Uribe Algorithm')
% disp(' A *E[xprime uprime] =B*[x u]')
% disp(' A = (dF/dxprimek dF/duprime), B =-(dF/dx dF/du)')
% A = [F1 F2]; B=[F3 F4]

% parameters which control numerical differentiation and Reiter solution
pnum     = p.NumWorkers;

packagesize = ceil(numstates/(3*pnum));
blocks      = ceil(numstates/packagesize);

% Absolute deviations
par.scaleval1 = 1e-5; %vector of numerical differentiation step sizes
par.scaleval2 = 1e-5; %vector of numerical differentiation step sizes

% jacobian wrt X', X

% disp('Computing Jacobian F1=DF/DXprime F3 =DF/DX')
% disp(['Total number of parallel blocks: ' num2str(blocks) '.'])
parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numstates);
    DF1 = zeros(length(Fb),length(range));
    DF3 = zeros(length(Fb),length(range));
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Xct = range
        X = zeros(numstates,1);
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

%
for i=1:ceil(numstates/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numstates);
    F1(:,range)=FF1{i};
    F3(:,range)=FF3{i};
end


% jacobian wrt Y'
packagesize = ceil(numcontrols/(3*pnum));
blocks=ceil(numcontrols/packagesize);
% disp('Computing Jacobian F2 - DF/DYprime')
% disp(['Total number of parallel blocks: ' num2str(blocks) '.'])
parfor bl=1:blocks
    range = ((bl-1)*packagesize+1):min(packagesize*bl,numcontrols);
    DF2=zeros(length(Fb),length(range));
    cc = zeros(numcontrols,1);
    ss = zeros(numstates,1);
    for Yct = range
        Y  = zeros(numcontrols,1);
        h = par.scaleval2;
        Y(Yct) = h;
        Fx = F(ss,ss,Y,cc);%#ok Fsys(NminusS_sp,X,C_sp,Cm_sp,Xss,Yss,Gamma_state,Gamma_control,InvGamma,par,mpar,grid,targets);
        DF2(:,Yct-(bl-1)*packagesize)=(Fx-Fb)/h;
    end
    FF{bl}=DF2;
%     disp(['Block number: ' num2str(bl) ' done.'])
end
%
for i=1:ceil(numcontrols/packagesize)
    range = ((i-1)*packagesize+1):min(packagesize*i,numcontrols);
    F2(:,range)=FF{i};
end

clear FF FF1 FF3;
% Derivative wrt Y (value functions today change only themselves)
%
cc = zeros(numcontrols,1);
ss = zeros(numstates,1);
for Yct = 0:oc-1
    Y  = zeros(numcontrols,1);
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
if nk>numstates
    if mpar.overrideEigen
        warning(['The Equilibrium is Locally Indeterminate, critical eigenvalue shifted to: ' num2str(ll(end-numstates))])
        slt = relev>ll(end-numstates);
        nk  = sum(slt);
    else
        error(['No Local Equilibrium Exists, last eigenvalue: ' num2str(ll(end-numstates))])
        %return
    end
elseif nk<numstates
    if mpar.overrideEigen
        warning(['No Local Equilibrium Exists, critical eigenvalue shifted to: ' num2str(ll(end-numstates))])
        slt = relev>ll(end-numstates);
        nk  = sum(slt);
    else
        error(['No Local Equilibrium Exists, last eigenvalue: ' num2str(ll(end-numstates))])
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
