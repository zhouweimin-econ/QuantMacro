% build transition matrix from states today to tomorrow 
% linear interpolation combining endogenous choices and exogenous shocks 

% Input variables:
%       anext_i   - nIx1 vector of chosen assets: anext_i(anow_i,idshock_i)
%       TransExo  - transition matrix for the exogenous shock     
%       nodesEndo - Grid points of the Endogenous Variables (e.g. Assets)
%       nS        - number of realziations of exogenous variables

function Pmat = BuildTrans(anext_i,TransExo,nodesEndo,nS)

% recover dimension
nI=length(nodesEndo);
amin=min(nodesEndo);
amax=max(nodesEndo);
nGrid=nS*nI;

% enlarge the choice vector for all possible future realizations of the shock
a_help=kron(anext_i,ones(nS,1));
a_help=min(max(a_help,amin),amax); % make sure we stay within the bounds of the Grid

% linear interpolation of a_help on the Grid
posEndo=lookupfor(nodesEndo,a_help,3);

fracEndo=1-(a_help-nodesEndo(posEndo))./(nodesEndo(posEndo+1)-nodesEndo(posEndo)); % interpolation weight (on left point)

%% weigths for possible realizations of the shock
WeightShock=reshape(kron(ones(nI,1),TransExo'),nI*nS^2,1);

rows=kron((1:nGrid)',ones(nS,1));
columns=kron(ones(nGrid,1),(0:(nS-1))'*nI);

P1=sparse(rows,columns+posEndo,fracEndo.*WeightShock,nGrid,nGrid);
P2=sparse(rows,columns+(posEndo+1),(1-fracEndo).*WeightShock,nGrid,nGrid);

Pmat=P1+P2;
