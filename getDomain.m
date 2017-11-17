function [ a, b, c, d, topoFunc ] = getDomain( testCase )

if strcmp( testCase, 'straka' ) || strcmp( testCase, 'strakaTopo' )
    a = -25600;
    b = 25600;
    c = 0;
    d = 6400;
    if strcmp( testCase, 'straka' )
        topoFunc = @(x)  zeros(size(x));
    else
%         topoFunc = @(x)  1000 * exp( -(32*(x-3200)/(b-a)).^2 );
        topoFunc = @(x)  500 * ( 1 - sin(2*pi*x./12800) );
    end
elseif strcmp( testCase, 'movingStraka' )
    a = -18000;
    b = 18000;
    c = 0;
    d = 6400;
%     topoFunc = @(x)  000 * exp( -(32*(x+000)/(b-a)).^2 );
    topoFunc = @(x)  500 * ( 1 - sin(2*pi*x./9000) );
elseif strcmp( testCase, 'doubleStraka' ) || strcmp( testCase, 'doubleStrakaTopo' ) ...
        || strcmp( testCase, 'doubleStrakaSmooth' ) || strcmp( testCase, 'doubleStrakaTopoSmooth' )
    a = -6400;
    b = 6400;
    c = 0;
    d = 6400;
    if strcmp( testCase, 'doubleStraka' ) || strcmp( testCase, 'doubleStrakaSmooth' )
        topoFunc = @(x)  zeros(size(x));
    elseif strcmp( testCase, 'doubleStrakaTopo' ) || strcmp( testCase, 'doubleStrakaTopoSmooth' )
        topoFunc = @(x)  1000 * exp( -(16*(x-1000)/(b-a)).^2 );
    end
elseif strcmp( testCase, 'bubble' ) || strcmp( testCase, 'bubbleSmooth' )
    a = 0;
    b = 10000;
    c = 0;
    d = 10000;
%     topoFunc = @(x)  000 * exp( -(8*(x-6000)/(b-a)).^2 );
    topoFunc = @(x)  500 * ( 1 + sin(2*pi*x./5000) );
elseif strcmp( testCase, 'mountainWaves' )
    a = 0;
    b = 10000;
    c = 0;
    d = 10000;
    topoFunc = @(x)  500 * exp( -(16*(x-5000)/(b-a)).^2 );
elseif strcmp( testCase, 'igw' )
    a = 0;
    b = 300000;
    c = 0;
    d = 10000;
    topoFunc = @(x)  zeros(size(x));
elseif strcmp( testCase, 'schar' )
    a = -25000;
    b = 25000;
    c = 0;
    d = 21000;
    hC = 250;
    lamC = 4000;
    aC = 5000;
    topoFunc = @(x)  hC * exp( -(x./aC) .^2 ) .* cos( pi*x./lamC ).^2;
end