function z = mq_x( d, x, y, rbforder )

z = x ./ sqrt( x.^2 + y.^2 + (rbforder*d).^2 );

z = z ./ d;