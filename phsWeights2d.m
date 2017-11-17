function [ W, W1, W2, Whv ] = phsWeights2d( X, Y, xe, ye, e1, e2, ...
	rbforder, mn, np, K, rbfType )

[ phi, phi_x, phi_y, phiHV ] = getFunctions( rbfType );

N = size(X,1);
n = size(X,2);

X = X - repmat(xe,1,n);
Y = Y - repmat(ye,1,n);
d = sqrt( X(:,n).^2 + Y(:,n).^2 );

X1 = X.*repmat(e1(:,1),1,n) + Y.*repmat(e1(:,2),1,n);
X2 = X.*repmat(e2(:,1),1,n) + Y.*repmat(e2(:,2),1,n);

W   = zeros( N, n );
W1  = zeros( N, n );
W2  = zeros( N, n );
Whv = zeros( N, n );

A = zeros( n+np, n+np );

for i = 1 : N
    
	x = X1(i,:);
	y = X2(i,:);
    
	xx = meshgrid(x);
	yy = meshgrid(y);
    A(1:n,1:n) = phi( d(i), xx.'-xx, yy.'-yy, rbforder );
    for j = 1 : np
        tmp = x.^mn(j,1) .* y.^mn(j,2) ./ d(i)^sum(mn(j,:));
        A(1:n,n+j) = tmp;
        A(n+j,1:n) = tmp;
    end
    
	b = [ phi(d(i),0-x,0-y,rbforder), zeros(1,np) ];
	b1 = [ phi_x(d(i),0-x,0-y,rbforder), zeros(1,np) ];
	b2 = [ phi_y(d(i),0-x,0-y,rbforder), zeros(1,np) ];
    bhv = [ phiHV(d(i),0-x,0-y,rbforder,K), zeros(1,np) ];
    if np>=1
        b(n+1) = 1;
    end
    if np>=3
        b1(n+2) = 1/d(i);
        b2(n+3) = 1/d(i);
    end
	if K==2 && np>=15
        bhv(n+11:n+15) = [ 24, 0, 8, 0, 24 ] ./ d(i)^4;
    elseif K==1 && np>=6
       	bhv(n+4:n+6) = [ 2, 0, 2 ] ./ d(i)^2;
	end
	w = [ b; b1; b2; bhv ] / A;
	w = w(:,1:n);
	W(i,:)   = w(1,:);
	W1(i,:)  = w(2,:);
	W2(i,:)  = w(3,:);
	Whv(i,:) = w(4,:);
end
