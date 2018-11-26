function [JD,EqQuant,c_star,m_star]  = denHaan_Reiter(JDminus,Sminus,...
    Control_F,Control_GUESS,ControlSS,indexMUdct,par,mpar,grid,meshes,P,aggrshock,oc)


%% Initializations
mutil = @(c)(1./(c.^par.xi));

% Number of states, controls
Ny = length(indexMUdct);
NN   = mpar.nm*mpar.nh; % Number of points in the full grid

%% Indexes for LHS/RHS
% Indexes for Controls
mutil_cind = 1:Ny;
Qind  = Ny+1;
Yind  = Ny+2;
Wind  = Ny+3;
Rind  = Ny+4;
Nind  = Ny+5;
Kind  = Ny+6;

%% Control Variables (Change Value functions according to sparse polynomial)
Control      = Control_F;
Controlminus = Control_GUESS;

Control(end-oc+1:end)       = ControlSS(end-oc+1:end) + Control_F(end-oc+1:end,:);
Controlminus(end-oc+1:end)  = ControlSS(end-oc+1:end) + (Control_GUESS(end-oc+1:end,:));

% Controls
XX             = zeros(NN,1);
XX(indexMUdct) = Control(mutil_cind);
% mutil_c        = idct(reshape(XX,[mpar.nm,mpar.nh]));
mutil_c        = (reshape(XX,[mpar.nm,mpar.nh]));
mutil_c        = (mutil(mutil_c(:)+ControlSS((1:NN))));

K  = exp(Control(Kind ));
Q  = exp(Control(Qind ));
R  = exp(Control(Rind ));

%% States
Hminus  = meshes.h(:)'*JDminus(:); 
Kminus  = meshes.m(:)'*JDminus(:);

%% Initialize
DD = @(Kguess)PriceD(Kguess,Q,R,Kminus,Hminus,Sminus,mutil_c,mpar,par,grid,meshes,aggrshock,P,JDminus);
QuantGuess=[K];

options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10);

[EqQuant]=fsolve(@(p)DD(p),QuantGuess,options);

[~,JD,c_star,m_star]= DD(EqQuant);
% fprintf('%2.2f \n',(abs(EqQuant-QuantGuess)./QuantGuess)*100);

end

function [QuantDiff,JD_new,c_star,m_star]= PriceD(QuantGuess,Q,R,Kminus,Hminus,Sminus,mutil_c,mpar,par,grid,meshes,aggrshock,P,JDminus)

K=QuantGuess(1);

%% Aggregate Output
switch(aggrshock)
    case('TFP')
        TFP=exp(Sminus);
end

Nminus    =  (TFP*par.alpha*Kminus.^(1-par.alpha)).^(1/(1-par.alpha+par.gamma));

% Wage Rate
Wminus =  TFP *par.alpha.* (Kminus./(Nminus)).^(1-par.alpha);
% Return on Capital
Rminus =  TFP *(1-par.alpha).* ((Nminus)./Kminus).^(par.alpha)  - par.delta;
%Price of Capital
Qminus=((K/Kminus)-1)*par.phi+1;

% Wages net of leisure services
WW=par.gamma/(1+par.gamma)*(Nminus./Hminus)*Wminus*ones(mpar.nm,mpar.nh);

%% Incomes (grids)
inc.labor   = WW.*(meshes.h);
inc.money   = (Rminus+Qminus)*meshes.m;
inc.profits  = 1/2*par.phi*((K-Kminus).^2)./Kminus;

%% Update policies

Raux = (R+Q)/(Qminus); 
EVm = reshape(reshape(Raux(:).*mutil_c,[mpar.nm mpar.nh])*P',[mpar.nm, mpar.nh]);

[c_star,m_star] = EGM_policyupdate(EVm,Rminus,Qminus,inc,meshes,grid,par,mpar);

%% Update distribution
% find next smallest on-grid value for money choices
weight11  = zeros(mpar.nm, mpar.nh,mpar.nh);
weight12  = zeros(mpar.nm, mpar.nh,mpar.nh);

% Adjustment case
[Dist_m,idm] = genweight(m_star,grid.m);

idm=repmat(idm(:),[1 mpar.nh]);
idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nh));

index11 = sub2ind([mpar.nm mpar.nh],idm(:),idh(:));
index12 = sub2ind([mpar.nm mpar.nh],idm(:)+1,idh(:));

for hh=1:mpar.nh
    
    %Corresponding weights
    weight11_aux = (1-Dist_m(:,hh));
    weight12_aux =  (Dist_m(:,hh));
    
    % Dimensions (mxk,h',h)
    weight11(:,:,hh)=weight11_aux(:)*P(hh,:);
    weight12(:,:,hh)=weight12_aux(:)*P(hh,:);
end

weight11=permute(weight11,[1 3 2]);
weight12=permute(weight12,[1 3 2]);

rowindex=repmat(1:mpar.nm*mpar.nh,[1 2*mpar.nh]);

H=sparse(rowindex,[index11(:); index12(:)],...
    [weight11(:); weight12(:)],mpar.nm*mpar.nh,mpar.nm*mpar.nh); % mu'(h',k'), a without interest


JD_new=JDminus(:)'*H;
JD_new=reshape(JD_new,[mpar.nm mpar.nh]);
%% Market clearing
KNEW  = meshes.m(:)'*JD_new(:);

QuantDiff=[KNEW-K];

end
function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end  % function
