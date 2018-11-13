function[contribution]=decomp(zm,sz1,sz2,sz3,sz4,sz5)
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
if (exist('sz3')~=1), error('Provide vector sz3 of a third shock'); end;
if (exist('sz4')~=1), error('Provide vector sz4 of a fourth shock'); end;
%if (exist('sz5')~=1), error('Provide vector sz5 of a fifth shock'); end;
szum=sz1+sz2+sz3+sz4;
g=zm-szum;
zm=g+szum;
VAR=var(g)+var(sz1)+var(sz2)+var(sz3)+var(sz4);
cg1=cov(g,sz1);
cg2=cov(g,sz2);
cg3=cov(g,sz3);
cg4=cov(g,sz4);
%cg5=cov(g,sz5);
c12=cov(sz1,sz2);
c13=cov(sz1,sz3);
c14=cov(sz1,sz4);
%c15=cov(sz1,sz5);
c23=cov(sz2,sz3);
c24=cov(sz2,sz4);
%c25=cov(sz2,sz5);
c34=cov(sz3,sz4);
%c35=cov(sz3,sz5);
%c45=cov(sz4,sz5);
COV=2*cg1(2,1)+2*cg2(2,1)+2*cg3(2,1)+2*cg4(2,1)+2*c12(2,1)+2*c13(2,1)+2*c14(2,1)...
    +2*c23(2,1)+2*c24(2,1)+2*c34(2,1);
var_zm=VAR+COV;
w_sz1=(var(sz1)+cg1(2,1)+c12(2,1)+c13(2,1)+c14(2,1))/var_zm;
w_sz2=(var(sz2)+cg2(2,1)+c12(2,1)+c23(2,1)+c24(2,1))/var_zm;
w_sz3=(var(sz3)+cg3(2,1)+c13(2,1)+c23(2,1)+c34(2,1))/var_zm;
w_sz4=(var(sz4)+cg4(2,1)+c14(2,1)+c24(2,1)+c34(2,1))/var_zm;
W=w_sz1^2+w_sz2^2+w_sz3^2+w_sz4^2;
wklad_sz1=w_sz1^2/W;
wklad_sz2=w_sz2^2/W;
wklad_sz3=w_sz3^2/W;
wklad_sz4=w_sz4^2/W;
contribution=[wklad_sz1 wklad_sz2 wklad_sz3 wklad_sz4]