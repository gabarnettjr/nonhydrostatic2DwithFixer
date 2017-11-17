function z = mq_y( d, x, y, rbforder )

z = y ./ sqrt( x.^2 + y.^2 + (rbforder*d).^2 );

z = z ./ d;