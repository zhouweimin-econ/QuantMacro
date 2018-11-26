% Aiyagari1 : SIMPLE MODEL OF CONSUMPTION AND SAVING UNDER UNCERTAINTY WITH CAPITAL
% Choose one of these three options to be equal to 1
% DEATON=1;
% HUGGETT=0;
% AIYAGARI=0;

clear
close all

% Choose one of these three options
DEATON=0;
HUGGETT=0;
AIYAGARI=1;


ns=300;
smin=-2;
smax=10;

BETA=0.94;
R=1.064;

S=linspace(smin,smax,ns);

N=2000; % Agents for simulations
T=500; % Periods for simulations

%-----------------------------------------------------
% Assign transition matrix and values for the income process
%-----------------------------------------------------

nz=2;
Z = [ 0.90 1.10 ]; 
P = [ 0.99 0.01
      0.01 0.99 ];

  
%-----------------------------------------------------
% Initialization of value functions
%-----------------------------------------------------

v=ones(nz,ns);
c=ones(1,ns);

V = ones(nz,ns)/(1-BETA);


if AIYAGARI==1;
    ALPHA=0.33;
    DELTAK=0.10;
    A=1;
    CAPITAL_SS = (ALPHA/(R-1+DELTAK))^(1/(1-ALPHA)) ;
    WAGE = A*(1-ALPHA)*(ALPHA*A/(R-1+DELTAK))^(ALPHA/(1-ALPHA)) ;
    W = WAGE*Z ;
end


%-----------------------------------------------------
% Iterate on value function until convergence
%-----------------------------------------------------

diffV   = 1;
iter    = 1;

% Initialize choices and value
newV = V;
EV = V;
idecS = zeros(nz,ns);
SPRIME = S ;


while (iter <= 500) & (diffV > 1e-6)

    % Calculate expected future value
    EV(:,1:ns)=P*V(:,1:ns) ;

%     EV1=P(1,:)*V(:,1:ns) ;
%     EV2=P(2,:)*V(:,1:ns) ;
%     EVV=[ EV1
%           EV2 ];
     
    for iz = 1:nz
        for is = 1:ns

            c = max( 1e-200, Z(iz) - SPRIME + R*S(is) ) ;
            v = log(c) + BETA*EV(iz,:) ;
 
            % "c", "v", EV(iz,:) are vectors of size [1,ns] 
            % We look in each state for the SPRIME that maximizes v

            [newV(iz,is), idecS(iz,is)] = max ( v ) ;

        end
    end

    
    % Use howard improvement algorithm to speed up calculations
    aiyagari1_hwd

    diffV = max(abs((newV(:)-V(:))./newV(:)));
    iter  = iter + 1;
    V     = newV ;
    disp(diffV);
    
    % Make plots only for illustrative purposes
%     subplot(2,1,1)
%     plot(S,V(1,:),'r'); 
%     title('Value function in bad state')
%     drawnow
%     subplot(2,1,2)
%     plot(S,V(2,:),'b'); 
%     title('Value function in good state')
%     drawnow;
    

end



%-----------------------------------------------------
% Calculate Decision Rules
%-----------------------------------------------------

Sdec = S(idecS(:,:)) ;

plot(S,Sdec(1,:),'r'); hold on; 
plot(S,Sdec(2,:),'b'); hold on;
plot(S,S,'k')
legend('Savings, low Z','Savings, high Z','45 line')




%-----------------------------------------------------
% Simulate consumption paths
%-----------------------------------------------------


if DEATON==1
    N=1;
    T=100;
    aiyagari1_sim
    break
end


aiyagari1_sim3
disp('Total savings')
disp(SS(end))
disp('Total capital')
if AIYAGARI==1;  disp(CAPITAL_SS); end
if HUGGETT==1; disp(0); CAPITAL_SS=0; end
if SS(end)>CAPITAL_SS
    disp('There is too much saving, interest rate must fall')
else
    disp('There is too little saving, interest rate must rise')
end
