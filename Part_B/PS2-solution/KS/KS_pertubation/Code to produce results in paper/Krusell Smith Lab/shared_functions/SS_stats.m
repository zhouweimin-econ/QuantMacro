% Steady state statistics    
    %%    
    
    clear targets
    grid.K    = grid.m*(sum(joint_distr,2));
    targets.K=grid.m*(sum(joint_distr,2));
    targets.KY=targets.K/Output;
    targets.Y=Output;
    BCaux_M=sum(sum(joint_distr,2),3);
    targets.m_bc=BCaux_M(1,:);
    par.RB = 1+R_fc;
    
    par.W=W_fc(1);
    par.N=N;
    par.R=R_fc;
    
    