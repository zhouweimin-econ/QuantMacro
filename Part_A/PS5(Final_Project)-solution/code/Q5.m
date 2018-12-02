clear;clf;close all;
cd '~/Desktop/PS5/code'  % in order to save png
% This file solves uncertainty case with CRRA
% Weimin Zhou
% codes modification based on:
% <George Hall, Econ 303: Advanced Macroeconomics I>
% <smodel2.m>: an Aiyagari (1994) QJE -style model. 
% http://people.brandeis.edu/~ghall/econ303/
% 
% revised after Luis's 2nd Lecture
%% II.5.1 baseline parameters values

sigma  = 2.0;             % risk aversion  
rho    = 0.06;            % Impatience coeff
beta   = 1/(1+rho);       % discount factor

R0     = 0.04;            % Initial guess of r
F      = -inf;            % Borrowing constraint parameter 


%sigma_y= 0.1; 			  % Uncertainty
sigma_y= 0.5;

gamma  = 0.95; 
prob     = [(1+gamma)/2,(1-gamma)/2;
        (1-gamma)/2,(1+gamma)/2];
% note:
% no delta, wage since in our model, no firms.
% hence, hereafter, only need to pin down r
%% II.5.2  Aiyagari (1994) parameters values
sigma    = 3;                          
beta     = 0.96;           
         
N        = 7;            
          
R0       = 0.04;         
F        = -10.0;  

sigma_y  = 0.4; 

gamma    = 0.36; 
prob	 = [(1+gamma)/2,(1-gamma)/2;
            (1+gamma)/2,(1-gamma)/2];
%% Discretization of the State Space
%  Endowment: 
Y   = [1-sigma_y,1-sigma_y];
nba = size(Y,2);
Y_high = Y(1,2);            % High value for income
Y_low  = Y(1,1);            % Low value for income

%   form asset grid  
maxaap    = 8;                     % Maximum value of asset grid   
incaap    = -5;                    % Minimum value of asset grid 


nbk      = 20;
%nbk   = round(maxaap/incaap+1); 
%intl    = 0.5;                    % Size of asset grid increments
intl	  = (maxaap-incaap)/(nbk-1); 

nasset    = round((maxaap-incaap)/intl+1); % number of grid points
Agrid     = [ incaap:intl:maxaap ]';

% pin down R util sum(lambda*A) = 0 
iter    = 1;
maxiter = 500;
toler   = 1e-4;
step    = 0.005;
R       = R0;
flag    = 1;

while  (flag ~= 0) && (iter <= maxiter);
   %  if we have firms, then:
   %  calculate rental rate of capital and wage
   %
   %  wage = (1-alpha) * A * K^(alpha)   * N^(-alpha);
   %  rent = (alpha)   * A * K^(alpha-1) * N^(1-alpha);
   % 
   %  tabulate the utility function such that for zero or negative
   %  consumption utility remains a large negative number so that
   %  such values will never be chosen as utility maximizing      
   %
   util1=-10000*ones(nasset,nasset);  % utility when high y    
   util2=-10000*ones(nasset,nasset);  % utility when low y
   
   for i=1:nasset
         asset=(i-1)*intl + incaap;
         
         % cons1:
         for j=1:nasset
               assetp = (j-1)*intl + incaap;
               cons = (1+R)*asset +Y_high-assetp;
               if assetp >= F && cons > 0;  % update utility only for positive consumption
                  util1(j,i)=(cons)^(1-sigma)/(1-sigma);
               end;
         
            cons2 = (1+R)*asset +Y_low-assetp;
               if assetp >= F && cons2 > 0;
                 util2(j,i)=(cons2)^(1-sigma)/(1-sigma);
               end;
         end;
   end;
   %
   %  initialize some variables
   %
   v       = zeros(nasset,2);
   decis   = zeros(nasset,2);
   test1    = 10;
   test2    = 10;
   [rs,cs] = size(util1); 
   %
   %  iterate on Bellman's equation and get the decision 
   %  rules and the value function at the optimum         
   %
   while (test1 ~= 0) || (test2 > .1);
       for i=1:cs;
           r1(:,i)=util1(:,i)+beta*(prob(1,1)*v(:,1)+ prob(1,2)*v(:,2));
           r2(:,i)=util2(:,i)+beta*(prob(2,1)*v(:,1)+ prob(2,2)*v(:,2));
       end;

       [tv1,tdecis1]=max(r1);
       [tv2,tdecis2]=max(r2);
       tdecis=[tdecis1' tdecis2'];
       tv=[tv1' tv2'];

       test1=max(any(tdecis-decis));
       test2=max(max(abs(tv - v))');
       v=tv;
       decis=tdecis;
   end;
   decis=(decis-1)*intl + incaap;  % Rescale the policy function
   %
   %   form transition matrix
   %   trans is the transition matrix from state at t (row)
   %   to the state at t+1 (column) 
   %   The eigenvector associated with the unit eigenvalue
   %   of trans' is  the stationary distribution. 
   %  
   g2=sparse(cs,cs); % Remember cs= nasset
   g1=sparse(cs,cs);
   
   for i=1:cs
       g1(i,tdecis1(i))=1;
       g2(i,tdecis2(i))=1;
   end
   trans=[ prob(1,1)*g1,prob(1,2)*g1; prob(2,1)*g2,prob(2,2)*g2];
   trans=trans';
 
   probst = (1/(2*nasset))*ones(2*nasset,1);
   test = 1;
   while test > 10^(-8);
       probst1  = trans*probst;
       test   = max(abs(probst1-probst));
       probst   = probst1;
   end; 
   %
   %   vectorize the decision rule to be conformable with probst
   %   calculate new aggregate assets  meanA
   %
   aa = decis(:);
   meanA = probst'*aa; %  Computing the expected a(R)
   %
   %  calculate measure over (k,s) pairs
   %  lambda has same dimensions as decis
   %
   lambda    = zeros(cs,2);
   lambda(:) = probst;
   %
   %   calculate stationary distribution of a
   %
   %[v1,d1]=eig(prob');
   %[dmax,imax]=max(diag(d1));
   %probst1=v1(:,imax);
   %ss=sum(probst1);
   %probst1=probst1/ss;
   proba = sum(lambda');     %  Stationary distribution of Assets 
   proba = proba';
   %
   %   form metric and update r (market clearing condition)
   %   if have K to pin down, then update K too.
   if iter == 1;
      A = meanA; 
      if meanA > 0.0;
         step = -step;
      end;
   end;
   Aold = A;
   Anew = meanA;
   if sign(Aold) ~= sign(Anew)
     step = -.5*step;
   end;
   if abs(step) >= toler;
      R = R + step;  % find market clearing condition and pick the MC R
   else
      flag = 0; % clear 
   end;
   
   A = Anew;
   iter = iter+1;
end;

ap		= Agrid(tdecis);
% Recovering consumption policy function
for j=1:nba;
   c(:,j)	= Y(j)+(1+R)*Agrid-ap(:,j);
end

plambda = lambda(:);  %  Stationary distribution of Assets 

% consumption distribution   
congood   = Y_high + (1+R)*Agrid - Agrid(tdecis(:,1));
conbad    = Y_low  + (1+R)*Agrid - Agrid(tdecis(:,2));
consum    = [congood conbad];
mean_con  = sum(diag(lambda'*consum)); %Multiplies the weight (lambda) by the value of consumption
mean_con2 = sum(diag(lambda'*(consum.^2)));
varcon    = (mean_con2 - mean_con^2);
[pconsum,index_c ] = sort(consum(:));


figure
plot(pconsum,plambda(index_c));
title('Consumption Distribution');
xlabel('consumption');
ylabel('population among all agents (%)');
saveas(gcf,'consumdist.png')

% utility(wealth) distribution 
UTILITY       = (consum).^(1-sigma)./(1-sigma);
mean_utility  = sum(diag(lambda'*UTILITY));
[putil,index_u ] = sort(UTILITY(:));

distribu        = 2*mean(reshape(probst,2,nbk));

figure
plot(Agrid,distribu)
title('Wealth Distribution')
ylabel('population among all agents (%)')
xlabel('asset holding')
saveas(gcf,'wealthdist.png')

%
%     calculate income distribution
%
income          = [((R+1)*Agrid + Y_high),((R+1)*Agrid + Y_low)]; 
mean_income     = sum(diag(lambda'*income));
[pinc,index_i ] = sort(income(:));
plambda = lambda(:);

figure
plot(pinc,plambda(index_i));
title('Income Distribution');
xlabel('income');
ylabel('population among all agents (%)');
saveas(gcf,'incomedist.png')

%% compare with KMP, compute quantile for my simulation (Q1-Q5) 
income_dist = plambda(index_i);
x_i2	= [];
for i	= 1:size(pinc,1)
   x_i2	= [x_i2 repmat(pinc(i),1,round(1000*income_dist(i),2))];
end

x_i2		 = x_i2'+1000;

edges_i = prctile(x_i2,[20,40,60,80,100]); 
for i = 1:length(edges_i)
   prcshare_i(i) = sum(x_i2<edges_i(i))/size(x_i2,1);
end

   
% Wealth
%distribu
x_w2		= [];
for i=1:size(putil,1)
   x_w2=[x_w2 repmat(putil(i),1,1000*round(probst(i),2))];
end

x_w2		 =	x_w2'+10000;
edges_w = prctile(x_w2,[20,40,60,80,100]); 
for i = 1:length(edges_w)
   prcshare_w(i) = sum(x_w2<=edges_w(i))/size(x_w2,1);
end

% Consumption
consum_dist		= plambda(index_c);
x_c2		= [];
for i=1:size(pconsum,1)
   x_c2=[x_c2 repmat(pconsum(i),1,1000*round(consum_dist(i),2))];
end

x_c2		 =	x_c2'+10000;
edges_c = prctile(x_c2,[20,40,60,80,100]); 
for i = 1:length(edges_c)
   prcshare_c(i) = sum(x_c2<edges_c(i))/size(x_c2,1);
end

disp('prcshare income')
prcshare_i
disp('--------------')

disp('prcshare wealth')
prcshare_w
disp('--------------')

disp('prcshare consumption')
prcshare_c
disp('--------------')