function [ N,R_fc,W_fc,Profits_fc,WW,RR,RBRB,Y ] = factor_returns(meshes,grid,par,mpar)
%factor_returns
  
    mc =  par.mu - (par.beta * log(par.PI) - log(par.PI))/par.kappa;
        
%%  GHH Preferences    
    N    =  (par.tau.*par.alpha.*grid.K.^(1-par.alpha).*mc).^(1/(1-par.alpha+par.gamma));   
    W_fc = par.alpha .*mc.* (grid.K./N).^(1-par.alpha);
    
    % Before tax return on capital
    R_fc = (1-par.alpha) .*mc.* (N./grid.K).^(par.alpha)- par.delta;
 
    Y = (N).^(par.alpha).*grid.K.^(1-par.alpha);
    Profits_fc = (1-mc)*Y - Y.*(1/(1-par.mu))./par.kappa./2 .*log(par.PI).^2;
 
    NW=par.gamma/(1+par.gamma).*(N/par.H).*W_fc;
    WW=NW*ones(mpar.nm,mpar.nk,mpar.nh); %Wages
    WW(:,:,end)=Profits_fc*par.profitshare;
    RR = R_fc; %Rental rates
    RBRB = par.RB/par.PI + (meshes.m<0).*(par.borrwedge/par.PI);


end

