function DCT=mydct(V,dim)
    n=size(V);
    dims=length(n);
    order=1:dims;
    neworder=order;
    neworder(dim)=[];
    neworder=[dim neworder];
    DC=mydctmx(n(dim));
    V=permute(V,neworder);
    nx=size(V);
    DCT = DC*reshape(V,[nx(1) numel(V)/nx(1)]);
    DCT = reshape(DCT,nx);
    DCT = permute(DCT,[2:dim 1 dim+1:dims]);
end