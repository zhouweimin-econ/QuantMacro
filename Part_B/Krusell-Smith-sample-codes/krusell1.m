%--------------------------------------------------------------------------
% SIMPLE MODEL OF CONSUMPTION AND SAVING UNDER IDIOSYNCRATIC AND AGGREAGATE
% UNCERTAINTY WITH CAPITAL (KRUSELL-SMITH)
%--------------------------------------------------------------------------


ns=180;
nK=3;
smin=-3;
smax=10;
BETA=0.90;

S=linspace(smin,smax,ns);

%-----------------------------------------------------
% Assign transition matrix and values for the income process
%-----------------------------------------------------

Zvec = [ 0.80 1.20 ]; nz=2;
PZ = [ 0.90 0.10
       0.10 0.90 ];

Avec =  [ 0.95 1.05 ]; na=2;
PA = [ 0.90 0.10
       0.10 0.90 ];

% Lump-together idiosyncratic and aggregate shock in a single P

[ AME ZME ]=ndgrid(Avec,Zvec);
P = kron(PA,PZ);

AMEvec = vec(AME) ;
ZMEvec = vec(ZME) ;

Yvec = AMEvec.*ZMEvec ;

ny = na*nz ;

disp('The transition matrix P')
disp(P)
disp(' ')
disp('The values of shocks')
disp(Yvec)
pause

%-----------------------------------------------------
% Initialize of value functions
%-----------------------------------------------------
v=ones(nK,ny,ns);
c=ones(1,ns);

V = ones(nK,ny,ns)/(1-BETA);


ALPHA=0.33;
DELTAK=0.10;
A=1;
CAPITAL_SS = (ALPHA/(1/BETA-1+DELTAK))^(1/(1-ALPHA)) ;
WAGE_SS = A*(1-ALPHA)*(ALPHA*A/(1/BETA-1+DELTAK))^(ALPHA/(1-ALPHA)) ;





% Initial vector for aggregate K
Kvec=linspace(0.50*CAPITAL_SS,2.00*CAPITAL_SS,nK);

Kmat = repmat(Kvec,[ny 1])
Kmat = Kmat';
AMEmat = repmat(AMEvec',[ nK 1 ])

pause

% The regression coefficients
% figure(gcf+1)

% For first iteration, use the coefficients below. For successive iterations, comment out 
% the lines below where values for b0, b1 and b2 are assigned and launch
% krusell1.m again

b0=0.0732;
b1=0.9695;
b2=0.9455;

disp(' ')
disp('The coefficients on K`=b0+b1*log(A)+b2*K')
disp([ b0 b1 b2 ])
Kmatprime = b0+b1*log(AMEmat)+b2*Kmat ;
disp(' ')

pause
 


for iK=1:nK
    for iy=1:ny
        Rmat(iK,iy)=1-DELTAK+ALPHA*AMEmat(iK,iy)/Kmat(iK,iy)^(1-ALPHA);
        WAGEmat(iK,iy) = (1-ALPHA)*AMEmat(iK,iy)*Kmat(iK,iy)^ALPHA ;
    end
end




%-----------------------------------------------------
% Iterate on value function until convergence
%-----------------------------------------------------
diffV   = 1;
iter    = 1;

% Initialize choices and value
newV = V;
EV = V;
idecS = zeros(nK,ny,ns);
SPRIME = S ;
v=ones(nK,ny,ns);



while (iter <= 200) & (diffV > 1e-6)

    % Interpolate value function on points outside the grid
    for iy=1:ny
        Vinterp(:,iy,1:ns)=interp1(Kmat(:,iy),V(:,iy,:),Kmatprime(:,iy),'linear','extrap');
    end

    % Calculate expected future value
    Vinterp_swap = permute(Vinterp,[2 1 3]);
    EV_swap = permute(EV,[2 1 3 ]);
    for is=1:ns
        EV_swap(:,:,is)=P*Vinterp_swap(:,:,is) ;
    end
    EV=permute(EV_swap,[2 1 3]);


    for iK = 1:nK
        for iy = 1:ny
            for is = 1:ns

                c = max( 1e-20, ZMEvec(iy)*WAGEmat(iK,iy) - SPRIME + Rmat(iK,iy)*S(is) )  ;

                % "c" and "v" are 1,ns vectors
                % For each state, look for the SPRIME that maximizes v

                v(iK,iy,:) = log(c') + BETA*squeeze(EV(iK,iy,:)) ;

                [newV(iK,iy,is), idecS(iK,iy,is)] = max ( v(iK,iy,:) ) ;
                

            end
        end
    end

    krusell1_hwd

    diffV = max(abs((newV(:)-V(:))./newV(:)));
    iter  = iter + 1;
    V     = newV ;
    disp(diffV);

end


%-----------------------------------------------------
% Calculate Decision Rules
%-----------------------------------------------------

Sdec = S(idecS(:,:,:)) ;

% figure(gcf+1)
% plot(S,squeeze(Sdec(3,1,:)),'r'); hold on
% plot(S,squeeze(Sdec(3,4,:)),'b'); hold on
% plot(S,S,'g')



%-----------------------------------------------------
% Simulate the economy
%-----------------------------------------------------

krusell1_sim