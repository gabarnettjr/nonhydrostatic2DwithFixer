function z = phs( d, x, y, rbforder )

if mod( rbforder, 2 ) == 0
    r = sqrt( x.^2 + y.^2 );
    z = r.^rbforder .* log(r);
    z(r<=eps) = 0;
else
    z = ( x.^2 + y.^2 ) .^ (rbforder/2);
end

z = 1./d.^rbforder .* z;