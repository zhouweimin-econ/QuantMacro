% Cakeeating1
% SOLVING BY VALUE FUNCTION ITERATION THE 
% SIMPLEST MODEL OF CONSUMPTION AND SAVING WITHOUT UNCERTAINTY 


clear
close all


ns=500;
smin=0;
smax=3;

BETA=0.90;

R=1.00; % An interest rate that generates going down to the lowest saving
Z=0.0;

S=linspace(smin,smax,ns);


%-----------------------------------------------------
% Initialization of value functions
% Initialize choices and value
%-----------------------------------------------------

v=ones(1,ns);
c=ones(1,ns);

% A decent guess
% V = log((1-BETA)*ones(1,ns))/(1-BETA);

V = zeros(1,ns);

% A bad guess
% V = (ones(1,ns))/(1-BETA);

newV = V;
EV = V;
idecS = zeros(1,ns);
SPRIME = S ;


%-----------------------------------------------------
% Iterate on value function until convergence
%-----------------------------------------------------
diffV   = 1;
iter    = 1;



while (iter <= 500) & (diffV > 1e-6)

    % Calculate expected future value
    EV=V ;

    for is = 1:ns

            c = max(1e-100, Z - SPRIME + R*S(is) ) ;
                        
            v = log(c) + BETA*EV ;
 
            % "c", "v", EV are vectors of size [1,ns] 
            % We look in each state for the SPRIME that maximizes v

            [newV(1,is), idecS(1,is)] = max ( v ) ;

    end
    
    
    diffV = max(abs((newV(:)-V(:))));
    iter  = iter + 1;
    V     = newV ;
    disp(diffV);
    
    plot(S,V); 
    pause(0.01)

end



%-----------------------------------------------------
% Calculate Decision Rules
%-----------------------------------------------------

figure(gcf+1)
Sdec = S(idecS(:)) ;


plot(S,Sdec,'r'); hold on; 
plot(S,S,'k')
legend('Savings','45 line')


%-----------------------------------------------------
% Compute and plot consumption function
%-----------------------------------------------------

figure(gcf+1)
for is=1:ns
Cdec(is) = Z - Sdec(is) + R*S(is) ;
end
plot(S,Cdec,'r'); hold on; 
title('Consumption function')


%-----------------------------------------------------
% Simulate an asset accumulation path
%-----------------------------------------------------

simS(1)=max(S);

for t = 2:15
    
    simS(t) = interp1(S,Sdec,simS(t-1)) ;
    simC(t) = Z - simS(t) + R*simS(t-1) ;
    simCOH(t) = Z + R*simS(t-1) ;

end

figure(gcf+1)
subplot(3,1,1)
plot(simS(2:end))
title('Simulated savings over time')
subplot(3,1,2)
plot(simC(2:end))
title('Simulated consumption')
subplot(3,1,3)
plot(simCOH(2:end))
title('Simulated beginning of period savings')
