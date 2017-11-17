function [ xx, zz, xxc, zzc, dx, dz ] = getNodes( topoFunc, a, b, c, d, nx, nz )

dx = (b-a) / (nx-1);
dz = (d-c) / (nz-1);

xx = linspace( a, b, nx ).';
xx = repmat( xx, 1, nz );
xxc = linspace( a+dx/2, b-dx/2, nx-1 ).';
xxc = repmat( xxc, 1, nz-1 );

zz = zeros( nx, nz );
for i = 1 : nx
    zz(i,:) = linspace( c+topoFunc(xx(i,1)), d, nz ) .';
end
zzc = ( zz(1:end-1,:) + zz(2:end,:) ) ./ 2;
zzc = ( zzc(:,1:end-1) + zzc(:,2:end) ) ./ 2;