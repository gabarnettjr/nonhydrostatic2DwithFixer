function z = ga( d, x, y, ep )

ep = ep ./ d;
z = exp( -ep.^2 .* ( x.^2 + y.^2 ) );