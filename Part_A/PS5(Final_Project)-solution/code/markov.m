% this function is attributed to Katrin Rabitsch
% obtained from:
% https://sites.google.com/site/katrinrabitsch/teaching/quantmacro2012
function [chain,state]=markov(T,n,s0,V);
%function [chain,state]=markov(T,n,s0,V);
%  chain generates a simulation from a Markov chain of dimension
%  the size of T
%
%  T is transition matrix
%  n is number of periods to simulate
%  s0 is initial state
%  V is the quantity corresponding to each state
%  state is a matrix recording the number of the realized state at time t
%
%

[r c]=size(T);
if nargin == 1;
  V=[1:r];
  s0=1;
  n=100;
end;
if nargin == 2;
  V=[1:r];
  s0=1;
end;
if nargin == 3;
  V=[1:r];
end;
%
if r ~= c;
  disp('error using markov function');
  disp('transition matrix must be square');
  return;
end;
%
for k=1:r;
  if sum(T(k,:)) ~= 1;
   disp('error using markov function')
    disp(['row ',num2str(k),' does not sum to one']);
    disp(' it sums to :'); 
    disp([ sum(T(k,:)) ]); 
    disp(['normalizing row ',num2str(k),'']);
    T(k,:)=T(k,:)/sum(T(k,:));
  end;
end;
[v1 v2]=size(V);
if v1 ~= 1 |v2 ~=r
  disp('error using markov function');
  disp(['state value vector V must be 1 x ',num2str(r),''])
  if v2 == 1 &v2 == r;
    disp('transposing state valuation vector');
    V=V';
  else;
    return;
  end;  
end
if s0 < 1 |s0 > r;
  disp(['initial state ',num2str(s0),' is out of range']);
  disp(['initial state defaulting to 1']);
  s0=1;
end;
%
%T
%rand('uniform');
X=rand(n-1,1);
s=zeros(r,1);
s(s0)=1;
cum=T*triu(ones(size(T)));
%
for k=1:length(X);
  state(:,k)=s;
  ppi=[0 s'*cum];

  s=((X(k)<=ppi(2:r+1)).*(X(k)>ppi(1:r)))';
    ppi1=ppi(2:r+1);
    ppi2=ppi(1:r);
    xa=X(k);
%     disp([ xa ]);
%     disp([ ppi1 ]);
%     disp([ ppi2 ]);
%     disp([ s ]); 
end;
chain=V*state;

