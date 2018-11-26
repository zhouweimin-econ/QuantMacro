function [DCT,DC]=mydct(V,dim,DC)


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
    DCT = DC*reshape(V,[nx(1) numel(V)/nx(1)]);
    DCT = reshape(DCT,nx);
    DCT = permute(DCT,[2:dim 1 dim+1:dims]);
end