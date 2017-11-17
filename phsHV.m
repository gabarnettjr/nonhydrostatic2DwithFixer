function z = phsHV( d, x, y, rbforder, K )

if mod( rbforder, 2 ) == 0

    r = sqrt( x.^2 + y.^2 );
    
    if K == 1
        z = rbforder * r.^(rbforder-2) .* ( 2 + rbforder*log(r) );
    elseif K == 2
        z = r.^(rbforder - 4).*rbforder.*(rbforder - 2).*(4*rbforder + rbforder^2*log(r) - 2*rbforder*log(r) - 4);
    elseif K == 3
        z = r.^(rbforder - 6).*rbforder.*(rbforder^2 - 6*rbforder + 8).*(rbforder^3*log(r) - 6*rbforder^2*log(r) - 24*rbforder + 8*rbforder*log(r) + 6*rbforder^2 + 16);
    elseif K == 4
        z = r.^(rbforder - 8)*rbforder*(rbforder^3 - 12*rbforder^2 + 44*rbforder - 48).*(176*rbforder + 44*rbforder^2*log(r) - 12*rbforder^3*log(r) + rbforder^4*log(r) - 48*rbforder*log(r) - 72*rbforder^2 + 8*rbforder^3 - 96);
    elseif K == 5
        z = r.^(rbforder - 10)*rbforder*(rbforder^4 - 20*rbforder^3 + 140*rbforder^2 - 400*rbforder + 384).*(140*rbforder^3*log(r) - 400*rbforder^2*log(r) - 1600*rbforder - 20*rbforder^4*log(r) + rbforder^5*log(r) + 384*rbforder*log(r) + 840*rbforder^2 - 160*rbforder^3 + 10*rbforder^4 + 768);
    end
    
    z(r<=eps) = 0;
    
else

	if K == 1
        z = rbforder^2 * (x.^2+y.^2) .^ ( (rbforder-2)/2 );
    elseif K == 2
        z = (rbforder-2)^2*rbforder^2 * (x.^2+y.^2) .^ ( (rbforder-4)/2 );
    elseif K == 3
        z = (rbforder-4)^2*(rbforder-2)^2*rbforder^2 * (x.^2+y.^2) .^ ( (rbforder-6)/2 );
    elseif K == 4
        z = (rbforder-6)^2*(rbforder-4)^2*(rbforder-2)^2*rbforder^2 * (x.^2+y.^2) .^ ( (rbforder-8)/2 );
    elseif K == 5
        z = (rbforder-8)^2*(rbforder-6)^2*(rbforder-4)^2*(rbforder-2)^2*rbforder^2 * (x.^2+y.^2) .^ ( (rbforder-10)/2 );
    end
    
end

z = 1./d.^rbforder .* z;