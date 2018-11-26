
% Given the deciosn rule so far..
Sdec = S(idecS(:,:)) ;

% construct a grid for income assets nz,ns
[ Zm Sm ] = ndgrid(Z,S);

% consumption is a grid too of dimension nz,ns
Cd = max( 1e-20, Zm - Sdec + R.*Sm ) ;

% evaluate utility
U =  log(Cd) ;

% newV=EV;


% Iterate over the policy function for 50 periods or so...
% that is, assume the same policy is used forever
% use this as new expected utility

for iHOWARD = 1:50

    for iz=1:nz
        for is=1:ns
            newV(iz,is) = U(iz,is) + BETA * EV(iz,idecS(iz,is)) ;
        end
    end

    % Calculate expected future value
    for is=1:ns
        EV(:,:)=P*newV(:,:) ;
    end

end