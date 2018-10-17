% ====
% Weimin Zhou
% Chebyshev Approximation
function Tx=chebyshev(x,n)
X=x(:);
lx=size(X,1);

    if n<0;
        error('n should be a positive integer');
    end
    
switch n;

    case 0;
        Tx=[ones(lx,1)];
    case 1;
        Tx=[ones(lx,1) X];
    otherwise
        Tx=[ones(lx,1) X];
    for i=3:n+1;
        Tx=[Tx 2*X.*Tx(:,i-1)-Tx(:,i-2)];
    end
end