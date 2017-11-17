function z = iqHV( d, x, y, ep, K )

ep = ep ./ d;
r = sqrt( x.^2 + y.^2 );
if K == 1
	z = 4 * ep.^2 .* ( (ep.*r).^2 - 1 ) ./ ( (ep.*r).^2 + 1 ) .^ 3;
elseif K == 2
	z = 64 * ep.^4 .* ( (ep.*r).^4 - 4*(ep.*r).^2 + 1 ) ./ ( (ep.*r).^2 + 1 ) .^5;
elseif K == 3
	z = 2304 * ep.^6 .* ( (ep.*r).^6 - 9*(ep.*r).^4 + 9*(ep.*r).^2 - 1 ) ./ ( (ep.*r).^2 + 1 ) .^ 7;
end

% r = sqrt( x.^2 + y.^2 );
% R = d .* R;
% if K == 1
	% z = 4 * ( r.^2 - R.^2 ) ./ ( r.^2 + R.^2 ) .^ 3;
% elseif K == 2
	% z = 64 * ( r.^4 - 4*r.^2.*R.^2 + R.^4 ) ./ ( r.^2 + R.^2 ) .^ 5;
% elseif K == 3
	% z = 2304 * ( r.^6 - 9*r.^4.*R.^2 + 9*r.^2.*R.^4 - R.^6 ) ./ ( r.^2 + R.^2 ) .^ 7;
% elseif K == 4
	% z = 147456 * ( r.^8 - 16*r.^6.*R.^2 + 36*r.^4.*R.^4 - 16*r.^2.*R.^6 + R.^8 ) ./ ( r.^2 + R.^2 ) .^ 9;
% end
% z = d.^2 .* z;