function [grid,P,bounds] = Tauchen(rho,N,sigma, mue, type)
%TAUCHEN Generates a discrete approximation to an AR-1 Process, following 
%   Tauchen (1987).
%   [grid,P] = Tauchen(rho,N,sigma, mue, type) 
%
%   returns a state vector GRID and a transition probability matrix P. 
%   
%   Input arguments:
%   rho:    autocorrelation coefficient
%   N:      number of gridpoints
%   sigma:  long-run variance
%   mue:    mean of the ar-1 process
%   type: {importance | equi | simple} - type of grid-Transition generating
%   algorithm.
%
%   importance: importance sampling. Each bin has probability 1/N to
%               realize
%   equi:       bin-centers are equi-spaced between +-3 std
%   simple:     like equi + Transition Probabilities are calculated without
%               using integrals.
%   simple importance: like simple but with grid from importance
%

%   Author: Christian Bayer, Uni Bonn, 03.05.2010

if nargin==2
    sigma=1;
    mue=0;
end
if nargin<=4
    type='importance';
end
aux=max(strcmp(type,{'importance','equi','simple','simple importance'}));
if aux(1)==0    
    type='importance';
    warning('TAUCHEN:NoOpt','No valid type set. Importance sampling used instead')
end
switch type
    case{'importance'} % Importance Sampling
                
        grid_probs=linspace(0,1,N+1);   % generate equi-likely bins

        bounds=norminv(grid_probs(1:end));     % corresponding bin bounds
        
        % replace (-)Inf bounds, by finite numbers
        bounds(1)=bounds(2)-99999999;
        bounds(end)=bounds(end-1)+99999999;
        
        % Calculate grid - centers
        grid = 1*N*( normpdf(bounds(1:end-1)) - normpdf(bounds(2:end)));
                    
        sigma_e=sqrt(1-rho^2);% Calculate short run variance
        P=NaN(N); % Initialize Transition Probability Matrix
        
        for i=1:floor((N-1)/2)+1 % Exploit Symmetrie to save running time
            for j=1:N
                P(i,j)=N*quadl(@(x)pr_ij(x,bounds(j),bounds(j+1),N,rho,sigma_e)...
                    ,bounds(i),bounds(i+1)); % Evaluate Integral
            end
        end
       % Exploit Symmetrie Part II
        P(floor((N-1)/2)+2:N,:)=P((ceil((N-1)/2):-1:1),end:-1:1);
        
  
   case{'equi'} % use +-3 std equi-spaced grid                      
        %% Equi-spaced
        step=6/(N-1);
        grid=-3:step:3;

        bounds=[-9999999 grid(1:end-1)+step/2 9999999]; % Bounds of Bins
      
        sigma_e=sqrt(1-rho^2);% Calculate short run variance
        P=NaN(N); % Initialize Transition Probability Matrix
       
        for i=1:floor((N-1)/2)+1 % Exploit Symmetrie to save running time
            for j=1:N
                P(i,j)=quadl(@(x)pr_ij(x,bounds(j),bounds(j+1),N,rho,sigma_e)...
                    ,bounds(i),bounds(i+1));% Evaluate integral
                
                % Normalize by Probability of the respective bin
                P(i,j)=P(i,j)/(normcdf(bounds(i+1))-normcdf(bounds(i)));
            end
        end
        % Exploit Symmetrie Part II
        P(floor((N-1)/2)+2:N,:)=P((ceil((N-1)/2):-1:1),end:-1:1); 
        
        % Make sure P is a Probability Matrix
 
        
    case{'simple'} % use simple Transition Probabilities
        % Generate Grid
        step=12/(N-1);
        grid = -6:step:6;
        bounds=[];
        sigma_e=sqrt(1-rho^2);% Calculate short run STD
        P=NaN(N); % Initialize Transition Probability Matrix
        for i=1:N
            P(i,1)=normcdf((grid(1)+step/2-rho*grid(i))/sigma_e);
            P(i,N)=1-normcdf((grid(N)-step/2-rho*grid(i))/sigma_e);
            for j=2:N-1
                P(i,j)=normcdf((grid(j)+step/2-rho*grid(i))/sigma_e) ...
                    -normcdf((grid(j)-step/2-rho*grid(i))/sigma_e);
            end
        end
      
    case{'simple importance'} % use simple Transition Probabilities
        % Generate Grid
        grid_probs=linspace(0,1,N+1);   % generate equi-likely bins
        bounds=norminv(grid_probs);     % corresponding bin bounds
        
        % Calculate grid - centers
        grid = N*( normpdf(bounds(1:end-1)) - normpdf(bounds(2:end)) );
        
        % replace (-)Inf bounds, by finite numbers
        bounds(1)=bounds(2)-99999999999;
        bounds(end)=bounds(end-1)+99999999999;
        
        sigma_e=sqrt(1-rho^2);% Calculate short run variance
        P=NaN(N); % Initialize Transition Probability Matrix
        for i=1:floor((N-1)/2+1)
            P(i,1) = normcdf((bounds(2)-rho*grid(i))/sigma_e);
            P(i,N) = 1 - normcdf((bounds(end-1)-rho*grid(i))/sigma_e);
            for j=2:N-1
                P(i,j)= normcdf((bounds(j+1)-rho*grid(i))/sigma_e) ...
                    - normcdf((bounds(j)-rho*grid(i))/sigma_e);
            end
        end
        P(floor((N-1)/2+2):end,:)=P(ceil((N-1)/2):-1:1,end:-1:1);
        
end
% Make sure P is a Probability Matrix
ps=sum(P,2);
P=P./(repmat(ps,1,N));

grid=grid*sqrt(sigma) + mue;


function pij=pr_ij(x,bound1,bound2,N,rho,sigma_e) % Pointwise likelihood from x to end in [bound1, bound2]
pij=normpdf(x) .* (normcdf((bound2 - rho.*x)./sigma_e) ...
    - normcdf((bound1 - rho.*x)./sigma_e));

