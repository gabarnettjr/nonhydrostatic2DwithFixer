function z = gaHV( d, x, y, ep, K )

ep = ep ./ d;
r = sqrt( x.^2 + y.^2 );

if K == 1
    z = 4*ep^2 * exp(-(ep*r).^2) .* ( (ep*r).^2 - 1 );
elseif K == 2
    z = 16*ep^4 * exp(-(ep*r).^2) .* ( (ep*r).^4 - 4*(ep*r).^2 + 2 );
elseif K == 3
    z = 64*ep^6 * exp(-(ep*r).^2) .* ( (ep*r).^6 - 9*(ep*r).^4 + 18*(ep*r).^2 - 6 );
elseif K == 4
    z = 256*ep^8 * exp(-(ep*r).^2) .* ( (ep*r).^8 - 16*(ep*r).^6 + 72*(ep*r).^4 - 96*(ep*r).^2 + 24 );
end