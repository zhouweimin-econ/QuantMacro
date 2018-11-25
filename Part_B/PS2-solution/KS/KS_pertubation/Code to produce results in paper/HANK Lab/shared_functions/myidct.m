function iDCT=myidct(V,dim,DC)
    n=size(V);
    dims=length(n);
    order=1:dims;
    neworder=order;
    neworder(dim)=[];
    neworder=[dim neworder];
    if nargin==2
    DC=mydctmx(n(dim));
    end
    V=permute(V,neworder);
    nx=size(V);
    iDCT = DC'*reshape(V,[nx(1) numel(V)/nx(1)]);
    iDCT = reshape(iDCT,nx);
    iDCT = permute(iDCT,[2:dim 1 dim+1:dims]);
end