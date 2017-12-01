function F = plottingFromSavedFiles

clc

%main parameters:
testCase = 'bubble';
rbfType = 'phs';
useMassFixer = 0;
nx = 200+1;
nz = 200+1;
rbforder    = 5;                    %the rbf parameter
polyorder   = 3;                    %for derivative approximations
n           = 11;                   %single-layer stencil size
sLayers     = 45/11;                %number of stencil layers
K           = 2;                    %HV parameter
gamma       = -2^-3;                %other HV parameter
tPlot       = 0 : 20 : 1500;

[ a, b, c, d, topoFunc ] = getDomain( testCase );

[ xx, zz, xxc, zzc, dx, dz ] = getNodes( topoFunc, a, b, c, d, nx, nz );

[ xxc, zzc, Nx, Nz, NxTop, NzTop, Tx, Tz, alp, bet, NC, bigTx, bigTz ] = ...
    addGhostNodes( xx, zz, xxc, zzc, a, b, n );

[ U, Cp, Cv, Rd, g, thetaBar, piBar, mu ] = getInitialConditions( testCase, xxc, zzc );

saveName = [ testCase, '/', rbfType, '_', ...
    'mf', num2str(useMassFixer), '_', ...
    'r', num2str(rbforder), 'p', num2str(polyorder), ...
    'n', num2str(n*sLayers), 'k', num2str(K), ...
    'g', num2str(abs(log2(abs(gamma)))), 'mu0', ...
    'dx', num2str((b-a)/(nx-1)), 'dz', num2str((d-c)/(nz-1)) ]

ind = getIndexes( a, b, xxc, zzc, n );

loops = length(tPlot);
F(loops) = struct( 'cdata', [], 'colormap', [] );

for i = 1 : length(tPlot)

    load( ['./matFiles/',saveName,'/',num2str(tPlot(i)),'.mat'], 'U' )
    
    fprintf( 1, '\nt=%g\nminPi=%g, maxPi=%g\nminU=%g, maxU=%g\nminW=%g, maxW=%g\nminTh=%g, maxTh=%g\n', ...
                tPlot(i), min(U(ind.m,1)-piBar(ind.m)), max(U(ind.m,1)-piBar(ind.m)), ...
                min(U(ind.m,2)), max(U(ind.m,2)), ...
                min(U(ind.m,3)), max(U(ind.m,3)), ...
                min(U(ind.m,4)), max(U(ind.m,4)) );

    plotContours( 1, testCase, xxc, zzc, xx, U, piBar, thetaBar, a, b, c, d, topoFunc, n, ind, Rd, Cp, Cv )
    
    F(i) = getframe;
    
end