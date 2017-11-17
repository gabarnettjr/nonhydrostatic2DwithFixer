function [ Wx, Wz, Whv, Wl ] = getWeights( computeWeights, simpleNeighbors, IDX, rad, XC, ZC, ...
    rbforder, n, K, sLayers, nx, nz, NC, mn, np, gamma, h, mu, rbfType )

[ phi, phi_x, phi_y, phiHV ] = getFunctions( rbfType );

nc = size(IDX,1);
stencilSize = round( n*sLayers );

Wx = zeros( stencilSize, nc );
Wz = zeros( stencilSize, nc );
Whv = zeros( stencilSize, nc );
Wl = zeros( stencilSize, nc );

if computeWeights == 1

    ii = repmat( 1:nc, stencilSize, 1 );
    jj = IDX.';

    A = zeros( stencilSize+np, stencilSize+np );
    b = zeros( 2, stencilSize+np );
    P = zeros( stencilSize, np );
    condA = zeros( nc, 1 );

    for i = 1 : nc

        if simpleNeighbors == 1
            xn = XC(i,:) - XC(i,1);
            zn = ZC(i,:) - ZC(i,1);
        else
            if sLayers == 3
                xn = XC(i,:) - XC(i,n+1);
                zn = ZC(i,:) - ZC(i,n+1);
            elseif sLayers == 5
                if i <= nx-1
                    xn = XC(i,:) - XC(i,n+1);
                    zn = ZC(i,:) - ZC(i,n+1);
                elseif i > (nz-2)*(nx-1)
                    xn = XC(i,:) - XC(i,3*n+1);
                    zn = ZC(i,:) - ZC(i,3*n+1);
                else
                    xn = XC(i,:) - XC(i,2*n+1);
                    zn = ZC(i,:) - ZC(i,2*n+1);
                end
            elseif sLayers == 7
                if i <= nx-1
                    xn = XC(i,:) - XC(i,n+1);
                    zn = ZC(i,:) - ZC(i,n+1);
                elseif i > nx-1 && i <=2*(nx-1)
                    xn = XC(i,:) - XC(i,2*n+1);
                    zn = ZC(i,:) - ZC(i,2*n+1);
                elseif i > (nz-3)*(nx-1) && i <= (nz-2)*(nx-1)
                    xn = XC(i,:) - XC(i,4*n+1);
                    zn = ZC(i,:) - ZC(i,4*n+1);
                elseif i > (nz-2)*(nx-1)
                    xn = XC(i,:) - XC(i,5*n+1);
                    zn = ZC(i,:) - ZC(i,5*n+1);
                else
                    xn = XC(i,:) - XC(i,3*n+1);
                    zn = ZC(i,:) - ZC(i,3*n+1);
                end
            end
        end

        xx = meshgrid(xn);
        zz = meshgrid(zn);

        %polynomial matrix:
        for j = 1 : np
            P(:,j) = xn.^mn(j,1) .* zn.^mn(j,2) ./ rad(i)^sum(mn(j,:));
        end
        A( 1:stencilSize, stencilSize+1:stencilSize+np ) = P;
        A( stencilSize+1:stencilSize+np, 1:stencilSize ) = P.';

        %differentiation weights for x and z:
		A(1:stencilSize,1:stencilSize) = phi( rad(i), xx.'-xx, zz.'-zz, rbforder );
		b(1,:) = [ phi_x(rad(i),0-xn,0-zn,rbforder), zeros(1,np) ];
		b(2,:) = [ phi_y(rad(i),0-xn,0-zn,rbforder), zeros(1,np) ];
        if np>=3
            b(1,stencilSize+2) = 1/rad(i);
            b(2,stencilSize+3) = 1/rad(i);
        end
        
        condA(i) = cond(A);
        
        w = b / A;
        Wx(:,i) = w(1,1:stencilSize);
        Wz(:,i) = w(2,1:stencilSize);

        bhv = [ phiHV(rad(i),0-xn,0-zn,rbforder,K), zeros(1,np) ];
		if K==2 && np>=15
			bhv(stencilSize+11:stencilSize+15) = [ 24, 0, 8, 0, 24 ] ./ rad(i)^4;
		elseif K==1 && np>=6
			bhv(stencilSize+4:stencilSize+6) = [ 2, 0, 2 ] ./ rad(i)^2;
		end
        w = bhv / A;
        Whv(:,i) = w(1:stencilSize);

        %laplacian weights:
%         A(1:stencilSize,1:stencilSize) = phi( rad(i), xx.'-xx, zz.'-zz, 3 );
%         bl = [ phiHV(rad(i),0-xn,0-zn,3,1), [0,0,0,2/rad(i)^2,0,2/rad(i)^2], zeros(1,np-6) ];
%         w = bl / A;
%         Wl(:,i) = w(1:stencilSize);

    end

    [ min(condA), max(condA) ]

    Wx = sparse( ii, jj, Wx, nc, NC, stencilSize*nc );
    Wz = sparse( ii, jj, Wz, nc, NC, stencilSize*nc );
    Whv = sparse( ii, jj, Whv, nc, NC, stencilSize*nc );
    Wl = sparse( ii, jj, Wl, nc, NC, stencilSize*nc );
    
    Whv = gamma * spdiags(h.^(2*K-1),-0:0,length(h),length(h)) * Whv;
    Wl = mu * Wl;
    
end