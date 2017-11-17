clc,clear

a  = -25600;
b  = 25600;
c  = 0;
d  = 6400;
nx = 64+1;
nz = 16+1;
topoFunc = @(x)  1000 * exp( -(16*(x+10000)/(b-a)).^2 );

dx = (b-a) / (nx-1);

[ xx, zz, xxc, zzc ] = getNodes( topoFunc, a, b, c, d, nx, nz, dx );

[ xx, zz, xxc, zzc, xNormals, zNormals ] = addGhostNodes( xx, zz, xxc, zzc, dx );

plotComputationalDomain( xx, zz, xxc, zzc, nx, nz, dx );