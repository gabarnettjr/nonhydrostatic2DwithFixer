function z = iq_y( d, x, y, ep )

ep = ep ./ d;
z = -2 * ep.^2 .* y ./ ( 1 + ep.^2 .* ( x.^2 + y.^2 ) ) .^ 2;

% z = -2*y ./ ( x.^2 + y.^2 + (d.*R).^2 ) .^ 2;
% z = d.^2 .* z;