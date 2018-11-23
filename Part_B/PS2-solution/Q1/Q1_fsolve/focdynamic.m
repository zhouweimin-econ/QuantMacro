function  F = focdynamic(consum,ap,a,Gamma,varphi,e,gamma,r0,w0)
% intratemporal condition 
F =  (1+r0)*a - ap + ((consum.^(-gamma).*e.*w0)/Gamma).^(varphi)*e*w0  - consum ;
end