function [ U, Cp, Cv, Rd, g, thetaBar, piBar, mu ] = getInitialConditions( testCase, xxc, zzc )

Cp = 1004;
Cv = 717;
Rd = Cp - Cv;
g = 9.81;

if strcmp( testCase, 'straka' ) || strcmp( testCase, 'strakaTopo' ) ...
        || strcmp( testCase, 'movingStraka' ) || strcmp( testCase, 'movingStrakaTopo' )
    
    thetaBar = 300 * ones(size(xxc));
    piBar = 1 - g ./ Cp ./ thetaBar .* zzc;
    xc = 0;
    zc = 3000;
    xr = 4000;
    zr = 2000;
    rTilde = sqrt( ((xxc-xc)./xr).^2 + ((zzc-zc)./zr).^2 );
    Tprime0 = zeros( size(rTilde) );
    ind = rTilde <= 1;
    Tprime0(ind) = -15/2 * ( 1 + cos(pi*rTilde(ind)) );
    thetaPrime0 = Tprime0 ./ piBar;
    
    pi0 = piBar(:);
    if strcmp( testCase, 'movingStraka' ) || strcmp( testCase, 'movingStrakaTopo' )
        u0 = 20*ones( size(pi0) );
    else
        u0 = zeros( size(pi0) );
    end
    w0 = zeros( size(pi0) );
    theta0 = thetaBar(:) + thetaPrime0(:);
    
    U = [ pi0, u0, w0, theta0 ];
    
    mu = 75;
    
elseif strcmp( testCase, 'doubleStraka' ) || strcmp( testCase, 'doubleStrakaTopo' ) ...
        || strcmp( testCase, 'doubleStrakaSmooth' ) || strcmp( testCase, 'doubleStrakaTopoSmooth' )
    
    thetaBar = 300 * ones(size(xxc));
    piBar = 1 - g ./ Cp ./ thetaBar .* zzc;
    xc1 = -6400;
    xc2 = 6400;
    zc = 3000;
    xr = 4000;
    zr = 2000;
    rTilde1 = sqrt( ((xxc-xc1)./xr).^2 + ((zzc-zc)./zr).^2 );
    rTilde2 = sqrt( ((xxc-xc2)./xr).^2 + ((zzc-zc)./zr).^2 );
    Tprime0 = zeros( size(rTilde1) );
    ind1 = rTilde1 <= 1;
    ind2 = rTilde2 <= 1;
    if strcmp( testCase, 'doubleStrakaSmooth' ) || strcmp( testCase, 'doubleStrakaTopoSmooth' )
        k = 2;
        Tprime0 = -15 * exp( -(k*rTilde1).^2 );
        Tprime0 = Tprime0 - 15 * exp( -(k*rTilde2).^2 );
    else
        Tprime0(ind1) = -15/2 * ( 1 + cos(pi*rTilde1(ind1)) );
        Tprime0(ind2) = -15/2 * ( 1 + cos(pi*rTilde2(ind2)) );
    end
    thetaPrime0 = Tprime0 ./ piBar;
    
    pi0 = piBar(:);
    u0 = zeros( size(pi0) );
    w0 = zeros( size(pi0) );
    theta0 = thetaBar(:) + thetaPrime0(:);
    
    U = [ pi0, u0, w0, theta0 ];
    
    mu = 75;
    
elseif strcmp( testCase, 'bubbleTopo' ) || strcmp( testCase, 'bubbleTopoSmooth' ) ...
        || strcmp( testCase, 'bubble' ) || strcmp( testCase, 'bubbleSmooth' )
    
    thetaBar = 300*ones(size(xxc));
    piBar = 1 - g ./ Cp ./ thetaBar .* zzc;
    R = 1500;
    xc = 5000;
    zc = 3000;
    r = sqrt( (xxc-xc).^2 + (zzc-zc).^2 );
    if strcmp( testCase, 'bubbleTopoSmooth' ) || strcmp( testCase, 'bubbleSmooth' )
        k = 1/900;
        thetaPrime0 = 2 * exp(-(k*r).^2);
    elseif strcmp( testCase, 'bubbleTopo' ) || strcmp( testCase, 'bubble' )
        ind = r < R;
        thetaPrime0 = zeros( size(r) );
        thetaPrime0(ind) = 2 * ( 1 - r(ind)./R );
    end
    
    pi0 = piBar(:);
    u0 = zeros( size(pi0) );
    w0 = zeros( size(pi0) );
    theta0 = thetaBar(:) + thetaPrime0(:);
    
    U = [ pi0, u0, w0, theta0 ];
    
    mu = 10;
    
elseif strcmp( testCase, 'mountainWaves' )
    
    thetaBar = 300 * ones(size(xxc));
    piBar = 1 - g ./ Cp ./ thetaBar .* zzc;
    
    pi0 = piBar(:);
    u0 = 5*ones( size(pi0) );
    w0 = zeros( size(pi0) );
    theta0 = thetaBar(:);
    
    U = [ pi0, u0, w0, theta0 ];
    
    mu = 0;
    
elseif strcmp( testCase, 'igw' ) || strcmp( testCase, 'igwTopo' )
    
    N = .01;
    theta0 = 300;
    thetaBar = theta0 * exp( (N^2/g) * zzc );
    
    piBar = (1-g^2./N^2./Cp./theta0) + g^2./N^2./Cp./theta0.*exp(-(N^2/g)*zzc);
    
    thetaC = .01;
    hC = 10000;
    aC = 5000;
    xC = 100000;
    thetaPrime0 = thetaC * sin( pi*zzc./hC ) ./ ( 1 + ((xxc-xC)./aC).^2 );
    
    pi0 = piBar(:);
    u0 = 20*ones(size(pi0));
    w0 = zeros(size(pi0));
    theta0 = thetaBar(:) + thetaPrime0(:);
    
    U = [ pi0, u0, w0, theta0 ];
    
    mu = 0;
    
elseif strcmp( testCase, 'schar' )
    
    N = .01;
    theta0 = 280;
    thetaBar = theta0 * exp( (N^2/g) * zzc );
    
    piBar = (1-g^2./N^2./Cp./theta0) + g^2./N^2./Cp./theta0.*exp(-(N^2/g)*zzc);
    
    pi0 = piBar(:);
    u0 = 10*ones(size(pi0));
    w0 = zeros(size(pi0));
    theta0 = thetaBar(:);
    
    U = [ pi0, u0, w0, theta0 ];
    
    mu = 0;
    
end

thetaBar = thetaBar(:);
piBar = piBar(:);
