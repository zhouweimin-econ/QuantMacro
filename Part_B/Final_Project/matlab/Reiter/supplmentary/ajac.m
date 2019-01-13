function y = ajac(a,d)
%AJAC Jacobian of vector field, either audi or function.
%
%   input audi array a:
%   y = AJAC(a)     audi matrix with Jacobian
%   y = AJAC(a,0) = aeval(ajac(a),0)
%
%   input function handle f:
%   g = AJAC(f,n) create a function handle g such that g(x1,...,xn) evaluates the 
%   Jacobian of f = f(x1,...,xn). Both input and output of g are of class double.
%   g = AJAC(f,var) same, but with an n-vector var of booleans. xi is treated as 
%   variable if var(i) = true and as parameter otherwise.
%
%   Example: <a href="matlab: ahelp(6)">Gaussian and mean curvature of a Klein bottle</a>
%   See also: acurl, adiv, agrad, ahess, alap
if isa(a,'function_handle')
  if isscalar(d)
    d = true(1,d);
  end
  [a1,a2,p] = prepfct(ainit([]),d,1);
  eval(['y = @(' a1 ') ajac(a(' a2 '),0);'])
  return
end
if nargin==2
  y = aeval(ajac(a),d);
  return
end
if size(a,2) > 1
  error('Evaluation for column vectors only.')
else
  y(numel(a),adim(a)) = audi();
  for i = 1:numel(a)
    y(i,:) = agrad(a(i))';
  end
end