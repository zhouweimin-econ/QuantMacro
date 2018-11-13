function[contribution]=decomp(zm,sz1,sz2)
%   put this routine in your current directory on MATLAB
%
%   out=decomp(zm,sz1,sz2,sz3,sz4,sz5)
%
%   Variance decomposition procedure.
%
%   Inputs:
%         zm - observable variable (vector maatrix) for which we want to decompose its variation
%         sz1, sz2, sz3, sz4, sz5 - shocks (vector matrices)
%         all input arguments should be of the same length
%   Outputs:
%         out - contribution of each shock's variation in the unconditional variability of the observable zm 
%
%   Created by K.Sobczak
%   January 2011
%
%   References:       Caselli, 2005; Growiec, 2010  
%
%   Please improve this routine in order to have any number of shocks !!! (with some loop I suppose)
if (exist('zm')~=1),  error('Provide vector zm of a variable for which you want to compute variance decomposition'); end;
if (exist('sz1')~=1), error('Provide vector sz1 of a first shock'); end;
if (exist('sz2')~=1), error('Provide vector sz2 of a second shock'); end;

%if (exist('sz5')~=1), error('Provide vector sz5 of a fifth shock'); end;
szum=sz1+sz2;
g=zm-szum;
zm=g+szum;
VAR=var(g)+var(sz1)+var(sz2);
cg1=cov(g,sz1);
cg2=cov(g,sz2);

c12=cov(sz1,sz2);

COV=2*cg1(2,1)+2*cg2(2,1)+2*c12(2,1);
var_zm=VAR+COV;
w_sz1=(var(sz1)+cg1(2,1)+c12(2,1))/var_zm;
w_sz2=(var(sz2)+cg2(2,1)+c12(2,1))/var_zm;

W=w_sz1^2+w_sz2^2;
wklad_sz1=w_sz1^2/W;
wklad_sz2=w_sz2^2/W;

contribution=[wklad_sz1 wklad_sz2]