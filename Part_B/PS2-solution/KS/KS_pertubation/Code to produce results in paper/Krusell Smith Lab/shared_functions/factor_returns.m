function [ N,W_fc,WW,RB,R_fc,Y ] = factor_returns(K,par,mpar)
%factor_returns
  
%%  GHH Preferences    
  
    N    =  (par.alpha.*K^(1-par.alpha)).^(1/(1-par.alpha+par.gamma));   
    W_fc =  par.alpha.* (K./N).^(1-par.alpha);        
    R_fc = (1-par.alpha).* (N./K).^(par.alpha)- par.delta;
    Y = N.^(par.alpha).*K.^(1-par.alpha);
 
    %%
    NW=par.gamma/(1+par.gamma).*(N/par.H).*W_fc;
    WW=NW*ones(mpar.nm,mpar.nh); %Wages
    RB = (1+R_fc);

end

