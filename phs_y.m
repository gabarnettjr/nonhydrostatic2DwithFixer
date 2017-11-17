function z = phs_y( d, x, y, rbforder )

if mod( rbforder, 2 ) == 0
    r = sqrt( x.^2 + y.^2 );
    z = y .* r.^(rbforder-2) .* ( 1 + rbforder*log(r) );
    z(r<=eps) = 0;
else
    z = rbforder .* y .* ( x.^2 + y.^2 ) .^ ((rbforder-2)/2);
end

z = 1./d.^rbforder .* z;