function J=jacobian(func,x0,param)
%
%  A function that uses numerical derivatives to compute 
%  the Jacobian matrix of any system of functions of several variables. 
%
% [Copyright: J. Ruiz]


auxi=diag(max(abs(x0)*1e-8,1e-8));
n=length(x0);
 for j=1:n
	J(:,j)=0.5.*(feval(func,x0+auxi(:,j),param)-...
        feval(func,x0-auxi(:,j),param))/auxi(j,j);
 end