function z = mq( d, x, y, rbforder )

z = sqrt( x.^2 + y.^2 + (rbforder.*d).^2 );

z = z ./ d;