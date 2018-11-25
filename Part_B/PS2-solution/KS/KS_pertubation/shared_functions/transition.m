function P=transition(N,rho,sigma_e,bounds)
% Calculate transition probability matrix for a given grid for a Markov
% chain with long-run variance equal to 1 and mean 0


P=NaN(N); % Initialize Transition Probability Matrix
for i=1:floor((N-1)/2)+1 % Exploit Symmetry to save running time
    for j=1:N
        P(i,j)=quadl(@(x)pr_ij(x,bounds(j),bounds(j+1),rho,sigma_e)...
            ,bounds(i),bounds(i+1));% Evaluate integral
        
        % Normalize by Probability of the respective bin
        P(i,j)=P(i,j)/(normcdf(bounds(i+1))-normcdf(bounds(i)));
    end
end
% Exploit Symmetry Part II
P(floor((N-1)/2)+2:N,:)=P((ceil((N-1)/2):-1:1),end:-1:1);
ps=sum(P,2);
P=P./(repmat(ps,1,N));

function pij=pr_ij(x,bound1,bound2,rho,sigma_e) % Pointwise likelihood from x to end in [bound1, bound2]
pij=normpdf(x) .* (normcdf((bound2 - rho.*x)./sigma_e) ...
    - normcdf((bound1 - rho.*x)./sigma_e));