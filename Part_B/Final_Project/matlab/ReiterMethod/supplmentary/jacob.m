% computes Jacobian by forward differences
% Input:
% func: string, name of function
% x: point at which to take Jacobian
% step: scalar, relative stepwidth;
% more arguments will be passed on to the function;
%
function jac = jacob(func,x,step,varargin)
  f0 = feval(func,x,varargin{:});
  n = size(x,1);
  m = size(f0,1);
  jac = zeros(m,n);
  x0 = x;
  for i=1:n
    step2 = step*max(1,x0(i));
    x = x0;
    x(i) = x0(i) + step2;
    jac(1:m,i) = (feval(func,x,varargin{:}) - f0)/step2;
  end;
