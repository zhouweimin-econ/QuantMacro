function [c_guess,k_star,joint_distr,R_fc,W_fc,Y,N,grid]=steadystate(P_H,grid,mpar,par,meshes)
  %     Prepare items for EGM
  %     Layout of matrices:
  %       Dimension 1: capital k
  %       Dimension 2: stochastic human capital h
          
  KL=0.5*grid.K;
  KH=2*grid.K;

  [excessL]=excessK(KL, grid,P_H,mpar,par,meshes);

  [excessH]=excessK(KH, grid,P_H,mpar,par,meshes);

  if sign(excessL)==sign(excessH)
    disp('ERROR! Sign not diff')
  end

  %% Brent
  fa=excessL; fb=excessH;
  a=KL;
  b=KH;

  if fa*fb>0
    error('f(a) und f(b) sollten unterschiedliche Vorzeichen haben');
  end

  c=a; fc=fa;   %Zu Beginn ist c = a

  c=a; fc=fa; d=b-a; e=d;

  iter=0;
  maxiter=1000;

  while iter<maxiter
    iter=iter+1;

    if fb*fc>0
      c=a; fc=fa; d=b-a; e=d;
    end

    if abs(fc)<abs(fb)
      a=b; b=c; c=a;
      fa=fb; fb=fc; fc=fa;
    end

    tol=2*eps*abs(b)+mpar.crit; m=(c-b)/2; %Toleranz

    if (abs(m)>tol) && (abs(fb)>0) %Verfahren muss noch durchgeführt werden

      if (abs(e)<tol) || (abs(fa)<=abs(fb))
        d=m; e=m;
      else
        s=fb/fa;
        if a==c
          p=2*m*s; q=1-s;
        else
          q=fa/fc; r=fb/fc;
          p=s*(2*m*q*(q-r)-(b-a)*(r-1));
          q=(q-1)*(r-1)*(s-1);
        end
        if p>0
          q=-q;
        else
          p=-p;
        end
        s=e; e=d;
        if ( 2*p<3*m*q-abs(tol*q) ) && (p<abs(s*q/2))
          d=p/q;
        else
          d=m; e=m;
        end
      end

      a=b; fa=fb;

      if abs(d)>tol
        b=b+d;
      else
        if m>0
          b=b+tol;
        else
          b=b-tol;
        end
      end
    else
      break;
    end

    [fb,c_guess]=excessK(b, grid,P_H,mpar,par,meshes);

  end

  Kcand=b;
  grid.K=b;

  %% Update

  [excess,c_guess,k_star,joint_distr,R_fc,W_fc,Y,N]=excessK(Kcand, grid,P_H,mpar,par,meshes);

end
function [excess,c_guess,k_star,joint_distr,R_fc,W_fc,Y,N]=excessK(K, grid,P_H,mpar,par,meshes)

    [ N,W_fc,WW,par.RB,R_fc,Y] = factor_returns(K,par,mpar);
    
    [c_guess,inc]=policyguess(meshes,WW,par);
    
    % Solve Policies and Joint Distribution
%     disp('Solving household problem by EGM')
    [c_guess,k_star,distPOL]=...
        policies_SS(c_guess, grid, inc,P_H,mpar,par);
        
%     disp('Calc Joint Distr')
    
    [joint_distr]=JDiteration(k_star,P_H,mpar,grid);
    joint_distr=reshape(joint_distr,[mpar.nm mpar.nh]);
    
    AggregateSavings=k_star(:)'*joint_distr(:);
    excess=AggregateSavings-K;
end
