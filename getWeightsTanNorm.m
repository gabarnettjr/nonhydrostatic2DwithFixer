function [ Wx, Wz, Whv, Wl ] = getWeightsTanNorm( computeWeights, simpleNeighbors, IDX, rad, XC, ZC, ...
    rbforder, n, K, sLayers, nx, nz, NC, mn, np, gamma, h, mu, e1, e2 )

% zS = topoFunc( xx(:,1) );
% zS = ( zS(1:end-1) + zS(2:end) ) / 2;
% zzcStar = ( zzc - repmat(zS,1,nz-1) ) ./ ( d - repmat(zS,1,nz-1) );
% figure(9),surf(xxc,zzc,zzcStar),axis([a,b,c,d,0,1]),pause

nc = size(IDX,1);
stencilSize = round( n*sLayers );

W1 = zeros( stencilSize, nc );
W2 = zeros( stencilSize, nc );
Whv = zeros( stencilSize, nc );
Wl = zeros( stencilSize, nc );

Wx = zeros( stencilSize, nc );
Wz = zeros( stencilSize, nc );

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

        tmpx = e1(i,1)*xn + e1(i,2)*zn;
        tmpz = e2(i,1)*xn + e2(i,2)*zn;
        xn = tmpx;
        zn = tmpz;

        xx = meshgrid(xn);
        zz = meshgrid(zn);

        %polynomial matrix:
        for j = 1 : np
            P(:,j) = xn.^mn(j,1) .* zn.^mn(j,2) ./ rad(i)^sum(mn(j,:));
        end
        A( 1:stencilSize, stencilSize+1:stencilSize+np ) = P;
        A( stencilSize+1:stencilSize+np, 1:stencilSize ) = P.';

        %differentiation weights for e1 and e2 (tangent and normal):
        A(1:stencilSize,1:stencilSize) = phi( rad(i), xx.'-xx, zz.'-zz, rbforder );
        condA(i) = cond(A);
        b(1,:) = [ phi_x(rad(i),0-xn,0-zn,rbforder), [0,1,0]./rad(i), zeros(1,np-3) ];
        b(2,:) = [ phi_y(rad(i),0-xn,0-zn,rbforder), [0,0,1]./rad(i), zeros(1,np-3) ];
        w = b / A;
        W1(:,i) = w(1,1:stencilSize);
        W2(:,i) = w(2,1:stencilSize);

        %hyperviscosity weights:
        A(1:stencilSize,1:stencilSize) = phi( rad(i), xx.'-xx, zz.'-zz, 2*K+1 );
        bhv = [ phiHV(rad(i),0-xn,0-zn,2*K+1,K), zeros(1,np) ];
        if K==2 && np>=15
            bhv(stencilSize+11:stencilSize+15) = [ 24, 0, 8, 0, 24 ] ./ rad(i)^4;
        end
        w = bhv / A;
        Whv(:,i) = w(1:stencilSize);

        %laplacian weights:
        A(1:stencilSize,1:stencilSize) = phi( rad(i), xx.'-xx, zz.'-zz, rbforder );
        bl = [ phiHV(rad(i),0-xn,0-zn,rbforder,1), [0,0,0,2,0,2]./rad(i)^2, zeros(1,np-6) ];
        w = bl / A;
        Wl(:,i) = w(1:stencilSize);

    end

    maxCondA = max(condA)

    W1 = sparse( ii, jj, W1, nc, NC, stencilSize*nc );
    W2 = sparse( ii, jj, W2, nc, NC, stencilSize*nc );
    Whv = sparse( ii, jj, Whv, nc, NC, stencilSize*nc );
    Wl = sparse( ii, jj, Wl, nc, NC, stencilSize*nc );

    Wx = spdiags(e1(:,1),-0:0,nc,nc)*W1 + spdiags(e2(:,1),-0:0,nc,nc)*W2;
    Wz = spdiags(e1(:,2),-0:0,nc,nc)*W1 + spdiags(e2(:,2),-0:0,nc,nc)*W2;

    Whv = gamma * spdiags(h.^(2*K-1),-0:0,length(h),length(h)) * Whv;
    Wl = mu * Wl;
    
end