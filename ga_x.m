function z = ga_x( d, x, y, ep )

ep = ep ./ d;
z = -2*ep^2*x .* exp( -ep.^2 .* ( x.^2 + y.^2 ) );