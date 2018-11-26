% Aiyagari_sim3.m: Simulate decision rules of agents in Deaton-Huggett-Aiyagari model
% enforcing Law of Large Numbers for each period t. Calls llnenforce.
% Speeding up simulations by calculating decision rules once and for all

% Initialization of series

simC    = ones(N,T) ;
simZ    = ones(N,T) ;
simS    = zeros(N,T) ;
sZ      = ones(N,T) ;

%--------------------------------------------------------------------
% With many agents, set half of them to be in the good state at time 1
% half of them to be in the bad state
%--------------------------------------------------------------------

if N>1
    sZ(:,1) = [ 1*ones(N/2,1)
                2*ones(N/2,1) ];
else
    sZ(:,1) = 1 ;
end

simZ(N,1) = Z(sZ(1)) ;

%--------------------------------------------------------------------
% Set initial savings to be equal to non-stochastic ss of K
%--------------------------------------------------------------------

if AIYAGARI==1
    simS(N,1) = CAPITAL_SS ;
else
    simS(N,1) = 0 ;
end




%--------------------------------------------------------------------
% Generate random numbers from uniform distribution to simulate Markov
% chain
%--------------------------------------------------------------------

rand('state',1); % reset the random number generator
rand_process=rand(N,T);

PLR = P^1000;

%--------------------------------------------------------------------
% For each agent, generate a sequence of income realizations
%--------------------------------------------------------------------
wait_time = waitbar(0,'Please wait, generating income realizations');


for i = 1:N

    for t = 2:T

        % Given a number from the uniform distribution...
        rp = rand_process(i,t) ;

        % Depending on the transition matrix, assign individual to a
        % particular state
        iz = 1;
        for xx=1:1:length(P)
            if rp>sum(P(sZ(i,t-1),1:xx)) ; iz = xx+1; end
        end

        sZ(i,t) = iz ;

    end

    llnenforce
    wb=waitbar(i/N,wait_time);

end

close(wb)


Sdec1 = Sdec(1,:) ;
Sdec2 = Sdec(2,:) ;

wait_time2 = waitbar(0,'Please wait, calculating decision rules');

for t = 2:T

        ind1=find(sZ(:,t)==1); % find index of agents who are in low income state
        simS(ind1,t) = interp1q(S',Sdec1',simS(ind1,t-1)) ;

        ind2=find(sZ(:,t)==2); % find index of agents who are in high income state
        simS(ind2,t) = interp1q(S',Sdec2',simS(ind2,t-1)) ;

        wb=waitbar(t/T,wait_time2);


end








%--------------------------------------------------------------------
% Once simulated income and savings are calculated
% generate consumption and a bunch of nice plots
%--------------------------------------------------------------------

simZ = Z(sZ) ;
simC(:,1) = simZ(:,1);
simC(:,2:T) = simZ(:,2:T) - simS(:,2:T) + R*simS(:,1:T-1);
SS=mean(simS);
SS_std=std(simS);

figure(gcf+1)
subplot(3,1,1) ; plot(simZ(1,:)); title('Simulated income, guy 1')
subplot(3,1,2) ; plot(simS(1,:)); title('Simulated savings, guy 1')
subplot(3,1,3) ; plot(simC(1,:)); title('Simulated consumption, guy 1')

if AIYAGARI==1 | HUGGETT==1
    figure(gcf+1)
    subplot(3,1,1) ; hist(simS(:,T),20); title('Time T distribution of savings')
    subplot(3,1,2) ; plot(SS); title('Average Savings over time')
    subplot(3,1,3) ; plot(SS_std); title('Standard Deviation of Savings over time')
end

close(wb)
