function [ errX, errZ ] = testSpatialDerivatives

clc

%main parameters:
useMQ = 0;
useMassFixer = 1;
testCase = 'strakaTopo';
nx = 8*16+1;
nz = 16+1;
rbforder    = 5;                    %the rbf parameter
polyorder   = 4;                    %for derivative approximations
n           = 9;                    %single-layer stencil size
sLayers     = 5;                    %number of stencil layers
K           = 2;                    %HV parameter
gamma       = -2^-3;                %other HV parameter
tPlot       = 0 : 10 : 1300;

simpleNeighbors = 1;
computeWeights = 1;
seeContours = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mn = getMN;  np = ( polyorder + 1 ) * ( polyorder + 2 ) / 2;

[ a, b, c, d, topoFunc ] = getDomain( testCase );

[ xx, zz, xxc, zzc, dx, dz ] = getNodes( topoFunc, a, b, c, d, nx, nz );

[ cellAreasExact, areaRatiosExact ] = getCellAreas( xx, zz );

[ xxc, zzc, Nx, Nz, NxTop, NzTop, Tx, Tz, alp, bet, NC, bigTx, bigTz ] = ...
    addGhostNodes( xx, zz, xxc, zzc, a, b, n );

tmpX = xxc( (n+1)/2:end-(n-1)/2, 2:end-1 );
tmpZ = zzc( (n+1)/2:end-(n-1)/2, 2:end-1 );
[ IDX, rad, XC, ZC, h ] = getNearestNeighbors( simpleNeighbors, xxc, zzc, tmpX, tmpZ, nx, nz, n, sLayers );

ind = getIndexes( a, b, xxc, zzc, n );

[ U, Cp, Cv, Rd, g, thetaBar, piBar, mu ] = getInitialConditions( testCase, xxc, zzc );

saveName = [ testCase, '/NEW_', ...
    'mq', num2str(useMQ), ...
    'mf', num2str(useMassFixer), ...
    'r', num2str(rbforder), 'p', num2str(polyorder), ...
    'n', num2str(n*sLayers), 'k', num2str(K), ...
    'g', num2str(abs(log2(abs(gamma)))), 'mu0', ...
    'dx', num2str((b-a)/(nx-1)), 'dz', num2str((d-c)/(nz-1)) ]

% [ e1, e2 ] = getNewCoordinates( xxc, zzc, n, dx );
% [ Wx, Wz, Whv, ~ ] = getWeightsTanNorm( computeWeights, simpleNeighbors, IDX, rad, XC, ZC, ...
    % rbforder, n, K, sLayers, nx, nz, NC, mn, np, gamma, h, mu, e1, e2 );
[ Wx, Wz, Whv, ~ ] = getWeights( computeWeights, simpleNeighbors, IDX, rad, XC, ZC, ...
    rbforder, n, K, sLayers, nx, nz, NC, mn, np, gamma, h, mu, useMQ );
	
	%extrapolation to ghost nodes (or bndry):
[ IDXb, radb, IDXt, radt, Xb, Zb, Xt, Zt, xb, zb, xt, zt ]  = getGhostNeighbors( xxc, zzc, n, sLayers, ind );
% if useMQ == 1
    % [ Wb, Wt, halfWb, halfWt ] = getExtrapolationWeights( radb, radt, Xb, Zb, Xt, Zt, xb, zb, xt, zt, ...
        % n, sLayers, rbforder, mn, 3, useMQ );
% else
    [ Wb, Wt, halfWb, halfWt ] = getExtrapolationWeights( radb, radt, Xb, Zb, Xt, Zt, xb, zb, xt, zt, ...
        n, sLayers, rbforder, mn, np, useMQ );
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check the error in approximating spatial derivatives:

% f = @(x,z)  cos(2*pi/6400*x) .* sin(2*pi/6400*z);
% f_x = @(x,z)  -2*pi/6400*sin(2*pi/6400*x) .* sin(2*pi/6400*z);
% f_z = @(x,z)  cos(2*pi/6400*x) * 2*pi/6400.*cos(2*pi/6400*z);

f = @(x,z)  exp( -((x+10000)/6400).^2 -((z-3200)/6400).^2 );
f_x = @(x,z)  -2*(x+10000)/6400/6400 .* f(x,z);
f_z = @(x,z)  -2*(z-3200)/6400/6400 .* f(x,z);

xtmp = xx(:,1);
ztmp = topoFunc(xtmp);

F = f(xxc(:),zzc(:));

F(ind.gb) = sum( Wb.*F(IDXb), 2 );
F(ind.gt) = sum( Wt.*F(IDXt), 2 );

app_x = Wx * F;
app_z = Wz * F;
app_hv = Whv * F;
xxc = xxc( (n+1)/2:end-(n-1)/2, 2:end-1 );
zzc = zzc( (n+1)/2:end-(n-1)/2, 2:end-1 );

if seeContours == 1
    
    figure(2),clf
    subplot(4,1,1)
        contourf( xxc, zzc, f(xxc,zzc), 10, 'lineStyle', 'none' )
        axis( 'equal', 'off' )
        colorbar
        colormap(parula(128))
        hold( 'on' )
        plot3( xtmp, ztmp, zeros(size(xtmp)), 'k' )
        plot3( [b,b,a,a], [c,d,d,c], 0*[1,1,1,1], 'k' )
        hold( 'off' )
        drawnow
    subplot(4,1,2)
        contourf( xxc, zzc, (reshape(app_x,size(xxc))-f_x(xxc,zzc))/max(max(abs(f_x(xxc,zzc)))), 10, 'lineStyle', 'none' )
        axis( 'equal', 'off' )
        colorbar
        hold( 'on' )
        plot3( xtmp, ztmp, zeros(size(xtmp)), 'k' )
        plot3( [b,b,a,a], [c,d,d,c], 0*[1,1,1,1], 'k' )
        hold( 'off' )
        drawnow
    subplot(4,1,3)
        contourf( xxc, zzc, (reshape(app_z,size(xxc))-f_z(xxc,zzc))/max(max(abs(f_z(xxc,zzc)))), 10, 'lineStyle', 'none' )
        axis( 'equal', 'off' )
        colorbar
        hold( 'on' )
        plot3( xtmp, ztmp, zeros(size(xtmp)), 'k' )
        plot3( [b,b,a,a], [c,d,d,c], 0*[1,1,1,1], 'k' )
        hold( 'off' )
        drawnow
    subplot(4,1,4)
        contourf( xxc, zzc, reshape(app_hv,size(xxc)), 10, 'lineStyle', 'none' )
        axis( 'equal', 'off' )
        colorbar
        hold( 'on' )
        plot3( xtmp, ztmp, zeros(size(xtmp)), 'k' )
        plot3( [b,b,a,a], [c,d,d,c], 0*[1,1,1,1], 'k' )
        hold( 'off' )
        drawnow
    
end

errX = norm( f_x(xxc(:),zzc(:)) - app_x, 2 ) ./ norm( f_x(xxc(:),zzc(:)), 2 );
errZ = norm( f_z(xxc(:),zzc(:)) - app_z, 2 ) ./ norm( f_z(xxc(:),zzc(:)), 2 );