
N=100;
T=100;
TDROP=20;
rand('state',1); % reset the random number generator


% Initialization of series
sA     = ones(1,T);
sZ(1:N,1) = ones(N,1);
sY = ones(N,T);

simS = ones(N,T)*CAPITAL_SS ;

simY = ones(N,T);

simZ = ones(N,T);

simK = ones(1,T)*CAPITAL_SS ;

a_rand = rand(1,T) ; % Aggregate shock is common to everybody
z_rand = rand(N,T) ;

t=1;
simZ(1:N,1) = Zvec(sZ(1:N,1)) ;
simY(:,1)=Yvec(sY(:,1));


wait_time = waitbar(0,'Please wait...');

for t = 2:T
    
    % Generate indices of aggregate shocks
    ra = a_rand(t) ;
    
    ia = 1;
    for xx=1:1:length(PA)
        if ra>sum(PA(sA(t-1),1:xx)) ; ia = xx+1; end
    end

    sA(t) = ia ;
        
    for i = 1:1:N

        % Generate indices of idiosyncratic shocks
        rz(i) = z_rand(i,t) ; % Random number for agent i at time t

        % Allocate random number to appropriate column of transition matrix
        iz = 1;
        for xx=1:1:length(PZ)
        if rz(i)>sum(PZ(sZ(i,t-1),1:xx)) ; iz = xx+1; end
        end

        sZ(i,t) = iz ;
        simZ(i,t) = Zvec(sZ(i,t));

        % From ia and iz get appropriate column of kronecker of PA*PZ for
        % aggregate-idio state mix
        iy = na*(iz-1)+ia ;
                   
        Sdec_squeezed = squeeze(Sdec(:,iy,:)) ;
               
        simS(i,t) = interpn(Kvec,S,Sdec_squeezed,simK(t-1),simS(i,t-1),'linear') ;
        simS(i,t) = max(min(simS(i,t),S(end)),S(1)) ;
        
    end

    
   
    % Calculate aggregate capital stock
    simK(t) = mean([ simS(:,t) ]) ;
    
    wb=waitbar(t/T,wait_time);

    
end



 
 
% Calculate simulated series
 
SSi     = simS(:,1:end-1);
SSinext = simS(:,2:end);

KKnext = simK(2:end);
KK     = simK(1:end-1);

ZZ     = simZ(:,2:end);

AA     = Avec(sA(2:end)); % Given timing, A begins in period 2. checked

RR      = ALPHA*AA./KK.^(1-ALPHA) + 1 - DELTAK ;
WW     = (1-ALPHA)*AA.*KK.^ALPHA ;


CCi = ZZ.*repmat(WW,[size(SSi,1) 1]) - SSinext + repmat(RR,[size(SSi,1) 1]).*SSi ;

CC = mean([ CCi ]) ;

pause

subplot(2,2,1)
plot(KKnext); title('K')
subplot(2,2,2)
plot(AA); title('A')
subplot(2,2,3)
plot(CC); title('C')
subplot(2,2,4)
plot(RR); title('R')

figure(gcf+1)
NGUYS=3;
for i=1:1:NGUYS
subplot(NGUYS,3,NGUYS*(i-1)+1)
plot(ZZ(i,:)); title('A'); title('Zi')
subplot(NGUYS,3,NGUYS*(i-1)+2)
plot(SSinext(i,:)); title('Si')
subplot(NGUYS,3,NGUYS*(i-1)+3)
plot(CCi(i,:)); title('Ci')
end
figure(gcf+1)


disp('Old coefficients, constant, K and log(A)')
disp([ b0 b1 b2 ])
disp(' ')


disp('Regression: K`=b0+b1*log(A)+b2*K')
YY=(KKnext');
XX=[ ones(numel(KKnext),1) log(AA') (KK')  ];
[BB,BINT,R2,RINT,STATS] = regress(YY(TDROP:end),XX(TDROP:end,:)) ;
 
 
disp('New coefficients')
disp(BB')
disp(' ')
disp('R2 of the regression')
disp(STATS(1))
disp(' ')
 
% Here we update the coefficients
 wnb=0.2;
 
b0=wnb*BB(1)+(1-wnb)*b0;
b1=wnb*BB(2)+(1-wnb)*b1;
b2=wnb*BB(3)+(1-wnb)*b2;

close(wb)
