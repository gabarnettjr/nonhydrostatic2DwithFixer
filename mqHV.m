function z = mqHV( d, x, y, rbforder, K )

if K == 1
    
    z = (2*(rbforder*d).^2 + x.^2 + y.^2)./((rbforder*d).^2 + x.^2 + y.^2).^(3/2);
    
elseif K == 2
    
    z = (- 8*(rbforder*d).^4 + 8*(rbforder*d).^2.*x.^2 + 8*(rbforder*d).^2.*y.^2 + x.^4 + 2*x.^2.*y.^2 + y.^4)./((rbforder*d).^2 + x.^2 + y.^2).^(7/2);
    
end

z = z ./ d;