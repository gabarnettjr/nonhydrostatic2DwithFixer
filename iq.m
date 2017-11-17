function z = iq( d, x, y, ep )

ep = ep ./ d;
z = 1 ./ ( 1 + ep.^2 .* ( x.^2 + y.^2 ) );

% z = 1 ./ ( x.^2 + y.^2 + (d.*R).^2 );
% z = d.^2 .* z;
