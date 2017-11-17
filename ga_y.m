function z = ga_y( d, x, y, ep )

ep = ep ./ d;
z = -2*ep^2*y .* exp( -ep.^2 .* ( x.^2 + y.^2 ) );