
nLayers = [ 4, 8, 16, 32, 64 ];
errX_c1 = zeros( size(nLayers) );
errZ_c1 = zeros( size(nLayers) );
errX_c2 = zeros( size(nLayers) );
errZ_c2 = zeros( size(nLayers) );

for i = 1 : length(nLayers)
    [ errX_c1(i), errZ_c1(i) ] = nonhydrostatic2D( nLayers(i), 2, 4, 5, 5 );
    [ errX_c2(i), errZ_c2(i) ] = nonhydrostatic2D( nLayers(i), 2, 4, 7, 5 );
end

lw = 1.5;
ms = 10;

figure(3),clf
plot( log10(nLayers), log10(errX_c1), '-o', 'lineWidth', lw, 'markerSize', ms )
hold( 'on' )
plot( log10(nLayers), log10(errZ_c1), '-o', 'lineWidth', lw, 'markerSize', ms )
plot( log10(nLayers), log10(errX_c2), '-o', 'lineWidth', lw, 'markerSize', ms )
plot( log10(nLayers), log10(errZ_c2), '-o', 'lineWidth', lw, 'markerSize', ms )
% plot( [1.4,1.8], [-4,-4-2*.4], '--' )
% plot( [1.4,1.8], [-6,-6-4*.4], '--' )
hold( 'off' )
ell = legend( 'xc1', 'zc1', 'xc2', 'zc2' );
drawnow