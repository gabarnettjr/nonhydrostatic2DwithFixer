function z = phs_x( d, x, y, rbforder )

if mod( rbforder, 2 ) == 0
    r = sqrt( x.^2 + y.^2 );
    z = x .* r.^(rbforder-2) .* ( 1 + rbforder*log(r) );
    z(r<=eps) = 0;
else
    z = rbforder .* x .* ( x.^2 + y.^2 ) .^ ((rbforder-2)/2);
end

z = 1./d.^rbforder .* z;