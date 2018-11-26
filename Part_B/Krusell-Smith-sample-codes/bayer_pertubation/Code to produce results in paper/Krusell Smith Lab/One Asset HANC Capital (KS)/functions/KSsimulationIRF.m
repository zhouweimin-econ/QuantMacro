function[K_sim, ss_sim]=KSsimulationIRF(k_star,P_H,joint_distr,grid,mesh,mpar,par)

    K_sim  = zeros(mpar.maxlag+mpar.T(2,1),1);
    ss_sim     = zeros(mpar.maxlag+mpar.T(2,1),1);
    K_sim(1)=grid.K(ceil(mpar.nK/2));
    grid.s=log(grid.s);
    
for t=1:mpar.maxlag+mpar.T(2,1)-1
    weight11  = zeros(mpar.nm, mpar.nh,mpar.nh);
    weight12  = zeros(mpar.nm, mpar.nh,mpar.nh);

        
    [~,idS]=histc(ss_sim(t),grid.s);
    idS(ss_sim(t)>=grid.s(end)) = mpar.ns-1;
    idS(ss_sim(t)<=grid.s(1))   = 1;
    
    distS=(ss_sim(t) - grid.s(idS))/(grid.s(idS+1) - grid.s(idS));
    
    k_star_aux= squeeze(k_star(:,idS,:,:)) + ...
        distS*(squeeze(k_star(:,idS+1,:,:))- squeeze(k_star(:,idS,:,:)));
        
    [~,idK]=histc(K_sim(t),grid.K);
    idK(K_sim(t)>=grid.K(end)) = mpar.nK-1;
    idK(K_sim(t)<=grid.K(1))   = 1;
    
    distK=(K_sim(t) - grid.K(idK))/(grid.K(idK+1) - grid.K(idK));
    
    k_star_aux= squeeze(k_star_aux(:,:,idK)) + ...
        distK*(squeeze(k_star_aux(:,:,idK+1))- squeeze(k_star_aux(:,:,idK)));
    
    
    [Dist_m,idm] = genweight(k_star_aux,grid.m);
    
    idm=repmat(idm(:),[1 mpar.nh]);
    idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nh));
    
    index11 = sub2ind([mpar.nm mpar.nh],idm(:),idh(:));
    index12 = sub2ind([mpar.nm mpar.nh],idm(:)+1,idh(:));
    
    for hh=1:mpar.nh
        
        %Corresponding weights
        weight11_aux = (1-Dist_m(:,hh));
        weight12_aux =  (Dist_m(:,hh));
        
        % Dimensions (mxk,h',h)
        weight11(:,:,hh)=weight11_aux(:)*P_H(hh,:);
        weight12(:,:,hh)=weight12_aux(:)*P_H(hh,:);
    end
    
    weight11=permute(weight11,[1 3 2]);
    weight12=permute(weight12,[1 3 2]);
    
    rowindex=repmat(1:mpar.nm*mpar.nh,[1 2*mpar.nh]);
    
    H=sparse(rowindex,[index11(:); index12(:)],...
        [weight11(:); weight12(:)],mpar.nm*mpar.nh,mpar.nm*mpar.nh); % mu'(h',k'), without interest
    
    joint_distr_next=joint_distr*H;
    joint_distr    = joint_distr_next;
    K_sim(t+1)=joint_distr*mesh.m(:);
    ss_sim(t+1) = ss_sim(t)*par.rhoS;
    if t==mpar.T(2,1)
        ss_sim(t+1)=par.sigmaS;
    end
        
    
end

end
function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end  % function