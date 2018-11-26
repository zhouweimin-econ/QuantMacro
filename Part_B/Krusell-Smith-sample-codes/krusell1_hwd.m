% Given the deciosn rule so far..
Sdec = S(idecS(:,:,:)) ;

Wm = repmat(WAGEmat,[1 1 ns]);
Rm = repmat(Rmat,[1 1 ns]);

Zm_temp = repmat(ZMEvec,[1 nK ns]);
Zm = permute(Zm_temp,[2 1 3]);

[ Km Ym Sm ] = ndgrid(Kvec,Yvec,S);

Cd = max( 1e-20, Wm.*Zm - Sdec + Rm.*Sm ) ;


U =  log(Cd) ;

newV=EV;
Vinterp=EV;

for iHOWARD = 1:50


    for iK=1:nK
        for iy=1:ny
            for is=1:ns
                newV(iK,iy,is) = U(iK,iy,is) + BETA * EV(iK,iy,idecS(iK,iy,is)) ;
            end
        end
    end



    for iy=1:ny
        Vinterp(:,iy,:)=interp1(Kmat(:,iy),newV(:,iy,:),Kmatprime(:,iy),'linear','extrap');
    end
    Vinterp_temp = permute(Vinterp,[2 1 3]);

    % Calculate expected future value
    for is=1:ns
        EV(:,:,is)=(P*Vinterp_temp(:,:,is))' ;
    end



end