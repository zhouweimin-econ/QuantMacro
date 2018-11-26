function [P_Rouw, z_Rouw] = rouwen(rho_Rouw, mu_uncond, sig_uncond, n_R)
%ROUWEN   Rouwenhorst's method (1995) to approximate an AR(1) process using 
%   a  finite state Markov process. 
%
%   For details, see Rouwenhorst, G., 1995: Asset pricing  implications of 
%   equilibrium business cycle models, in Thomas Cooley (ed.), Frontiers of 
%   Business Cycle Research, Princeton University Press, Princeton, NJ.
% 
%   Suppose we need to approximate the following AR(1) process:
%
%                   y'=rho_Rouw*y+e
%
%   where abs(rho_Rouw)<1, sig_uncond=std(e)/sqrt(1-rho_Rouw^2) and 
%   mu_uncond denotes E(y), the unconditional mean of y. Let n_R be the 
%   number of grid points. n_R must be a positive integer greater than one.  
%
%   [P_Rouw, z_Rouw] = rouwen(rho_Rouw, mu_uncond, sig_uncond, n_R) returns  
%   the discrete state space of n_R grid points for y, z_Rouw, and 
%   the centrosymmetric transition matrix P_Rouw. Note that
%
%       1. z_Rouw is a column vector of n_R real numbers. 
%       2. The (i,j)-th element of P_Rouw is the conditional probability 
%          Prob(y'=z_Rouw(i)|y=z_Rouw(j)), i.e.
%
%                 P_Rouw(i,j)=Prob(y'=z_Rouw(i)|y=z_Rouw(j))
%
%           where z_i is the i-th element of vector z_Rouw. Therefore 
%
%           P_Rouw(1,j)+P_Rouw(2,j)+ ... +P_Rouw(n,j)=1 for all j.
%   
%   See also HITM_Z and HITM_S on how to simulate a Markov processes using 
%   a transition matrix and the grids. 
%
%   Damba Lkhagvasuren, June 2005

% CHECK IF abs(rho)<=1 
if abs(rho_Rouw)>1
    error('The persistence parameter, rho, must be less than one in absolute value.');
end

% CHECK IF n_R IS AN INTEGER GREATER THAN ONE.
if n_R <1.50001 %| mod(n_R,1)~=0 
    error('For the method to work, the number of grid points (n_R) must be an integer greater than one.');  
end

% CHECK IF n_R IS AN INTEGER.
if mod(n_R,1)~=0 
    warning('the number of the grid points passed to ROUWEN is not an integer. The method rounded n_R to its nearest integer.')
    n_R=round(n_R);
    disp('n_R=');
    disp(n_R);  
end

% GRIDS
step_R = sig_uncond*sqrt(n_R - 1); 
z_Rouw=[-1:2/(n_R-1):1]';
z_Rouw=mu_uncond+step_R*z_Rouw;

% CONSTRUCTION OF THE TRANSITION PROBABILITY MATRIX
p=(rho_Rouw + 1)/2;
q=p;

P_Rouw=[ p  (1-p);
        (1-q) q];
    
    for i_R=2:n_R-1
    a1R=[P_Rouw zeros(i_R, 1); zeros(1, i_R+1)];
    a2R=[zeros(i_R, 1) P_Rouw; zeros(1, i_R+1)];
    a3R=[zeros(1,i_R+1); P_Rouw zeros(i_R,1)];
    a4R=[zeros(1,i_R+1); zeros(i_R,1) P_Rouw];
    P_Rouw=p*a1R+(1-p)*a2R+(1-q)*a3R+q*a4R;
    P_Rouw(2:i_R, :) = P_Rouw(2:i_R, :)/2;
    end
    
P_Rouw=P_Rouw';

for i_R = 1:n_R
    P_Rouw(:,i_R) = P_Rouw(:,i_R)/sum(P_Rouw(:,i_R));
end
P_Rouw=P_Rouw';