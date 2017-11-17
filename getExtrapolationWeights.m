function [ Wb, Wt, halfWb, halfWt ] = getExtrapolationWeights( radb, radt, Xb, Zb, Xt, Zt, xb, zb, xt, zt, ...
    n, sLayers, rbforder, mn, np, rbfType )
	
[ phi, phi_x, phi_y, phiHV ] = getFunctions( rbfType );

stencilSize = round( n * sLayers );

%weights for extrapolating to the ghost nodes:

Wb = zeros( size(Xb) );
for i = 1 : size(Xb,1)
    xn = Xb(i,:) - xb(i);
    zn = Zb(i,:) - zb(i);
    xx = meshgrid(xn);
    zz = meshgrid(zn);
    A(1:stencilSize,1:stencilSize) = phi( radb(i), xx.'-xx, zz.'-zz, rbforder );
    for j = 1 : np
        tmp = xn.^mn(j,1) .* zn.^mn(j,2) ./ radb(i)^sum(mn(j,:));
        A(1:stencilSize,stencilSize+j) = tmp;
        A(stencilSize+j,1:stencilSize) = tmp;
    end
    b = [ phi(radb(i),0-xn,0-zn,rbforder), [1,zeros(1,np-1)] ];
    w = b / A;
    Wb(i,:) = w(1:stencilSize);
end

Wt = zeros( size(Xt) );
for i = 1 : size(Xt,1)
    xn = Xt(i,:) - xt(i);
    zn = Zt(i,:) - zt(i);
    xx = meshgrid(xn);
    zz = meshgrid(zn);
    A(1:stencilSize,1:stencilSize) = phi( radt(i), xx.'-xx, zz.'-zz, rbforder );
    for j = 1 : np
        tmp = xn.^mn(j,1) .* zn.^mn(j,2) ./ radt(i)^sum(mn(j,:));
        A(1:stencilSize,stencilSize+j) = tmp;
        A(stencilSize+j,1:stencilSize) = tmp;
    end
    b = [ phi(radt(i),0-xn,0-zn,rbforder), [1,zeros(1,np-1)] ];
    w = b / A;
    Wt(i,:) = w(1:stencilSize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%weights for extrapolating to the boundary:

halfWb = zeros( size(Xb) );
for i = 1 : size(Xb,1)
    xn = Xb(i,:) - (xb(i)+Xb(i,1))/2;
    zn = Zb(i,:) - (zb(i)+Zb(i,1))/2;
    xx = meshgrid(xn);
    zz = meshgrid(zn);
    A(1:stencilSize,1:stencilSize) = phi( radb(i), xx.'-xx, zz.'-zz, rbforder );
    for j = 1 : np
        tmp = xn.^mn(j,1) .* zn.^mn(j,2) ./ radb(i)^sum(mn(j,:));
        A(1:stencilSize,stencilSize+j) = tmp;
        A(stencilSize+j,1:stencilSize) = tmp;
    end
    b = [ phi(radb(i),0-xn,0-zn,rbforder), [1,zeros(1,np-1)] ];
    w = b / A;
    halfWb(i,:) = w(1:stencilSize);
end

halfWt = zeros( size(Xt) );
for i = 1 : size(Xt,1)
    xn = Xt(i,:) - (xt(i)+Xt(i,1))/2;
    zn = Zt(i,:) - (zt(i)+Zt(i,1))/2;
    xx = meshgrid(xn);
    zz = meshgrid(zn);
    A(1:stencilSize,1:stencilSize) = phi( radt(i), xx.'-xx, zz.'-zz, rbforder );
    for j = 1 : np
        tmp = xn.^mn(j,1) .* zn.^mn(j,2) ./ radt(i)^sum(mn(j,:));
        A(1:stencilSize,stencilSize+j) = tmp;
        A(stencilSize+j,1:stencilSize) = tmp;
    end
    b = [ phi(radt(i),0-xn,0-zn,rbforder), [1,zeros(1,np-1)] ];
    w = b / A;
    halfWt(i,:) = w(1:stencilSize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%