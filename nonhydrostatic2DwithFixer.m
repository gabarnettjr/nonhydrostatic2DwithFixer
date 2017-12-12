function nonhydrostatic2DwithFixer

clc

%main parameters:
testCase = 'doubleStrakaTopoSmooth';
rbfType = 'phs';
useMassFixer = 0;
nx = 512+1;
nz = 256+1;
rbforder    = 5;                    %the rbf parameter
polyorder   = 3;                    %for derivative approximations
n           = 11;                   %single-layer stencil size
sLayers     = 100/11;               %number of stencil layers
K           = 2;                    %HV parameter
gamma       = -2^-3;                %other HV parameter
tPlot       = 0 : 10 : 900;

dt = 1/24;
rkStages = 3;
simpleNeighbors = 1;      %if 1, just get neighbors instead of using layers
adamsBashforth = 0;
rbfMultistep   = 0;

%other parameters:
seeDomain      = 0;
seeStencils    = 0;
computeWeights = 1;
timeStep       = 1;
seeContours    = 0;
saveResults    = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

t = tPlot(1) : dt : tPlot(end);

mn = getMN;  np = ( polyorder + 1 ) * ( polyorder + 2 ) / 2;

[ a, b, c, d, topoFunc ] = getDomain( testCase );

[ xx, zz, xxc, zzc, dx, dz ] = getNodes( topoFunc, a, b, c, d, nx, nz );

[ cellAreasExact, areaRatiosExact ] = getCellAreas( xx, zz );

[ xxc, zzc, Nx, Nz, NxTop, NzTop, Tx, Tz, TxTop, TzTop, ...
	alp, bet, NC, bigTx, bigTz, bigNx, bigNz ] = addGhostNodes( xx, zz, xxc, zzc, a, b, n );

tmpX = xxc( (n+1)/2:end-(n-1)/2, 2:end-1 );
tmpZ = zzc( (n+1)/2:end-(n-1)/2, 2:end-1 );
[ IDX, rad, XC, ZC, h ] = ...
    getNearestNeighbors( simpleNeighbors, xxc, zzc, tmpX, tmpZ, nx, nz, n, sLayers, dx, dz );

% areaRatios = h/sum(h);
% cellAreas = areaRatios .* sum(cellAreasExact);
areaRatios = areaRatiosExact;
cellAreas = cellAreasExact;

ind = getIndexes( a, b, xxc, zzc, n );

plotComputationalDomain( seeDomain, xx, zz, xxc, zzc, nx, nz, dx, n, ind );

plotStencils( simpleNeighbors, seeStencils, xxc, zzc, Nx, Nz, NxTop, NzTop, alp, bet, XC, ZC, n, sLayers, nx, nz );

[ U, Cp, Cv, Rd, g, thetaBar, piBar, mu ] = getInitialConditions( testCase, xxc, zzc );

saveName = [ testCase, '/', rbfType, '_', ...
    'mf', num2str(useMassFixer), '_', ...
    'r', num2str(rbforder), 'p', num2str(polyorder), ...
    'n', num2str(n*sLayers), 'k', num2str(K), ...
    'g', num2str(abs(log2(abs(gamma)))), 'mu0', ...
    'dx', num2str((b-a)/(nx-1)), 'dz', num2str((d-c)/(nz-1)) ]

% [ e1, e2 ] = getNewCoordinates( xxc, zzc, n, dx );
% [ Wx, Wz, Whv, ~ ] = getWeightsTanNorm( computeWeights, simpleNeighbors, IDX, rad, XC, ZC, ...
%     rbforder, n, K, sLayers, nx, nz, NC, mn, np, gamma, h, mu, e1, e2 );
[ Wx, Wz, Whv, ~ ] = getWeights( computeWeights, simpleNeighbors, IDX, rad, XC, ZC, ...
    rbforder, n, K, sLayers, nx, nz, NC, mn, np, gamma, h, mu, rbfType );

%extrapolation to ghost nodes (or bndry):
[ IDXb, IDXt, Xb, Zb, Xt, Zt, xb, zb, xt, zt ]  = ...
	getGhostNeighbors( xxc, zzc, n, sLayers, ind, dx, dz );
[ Wb, ~, ~, ~ ] = phsWeights2d( Xb, Zb, xb, zb, [Tx,Tz], [Nx,Nz], ...
	rbforder, mn, np, K, rbfType );
[ halfWb, ~, ~, ~ ] = phsWeights2d( Xb, Zb, (xb+Xb(:,1))/2, (zb+Zb(:,1))/2, [Tx,Tz], [Nx,Nz], ...
	rbforder, mn, np, K, rbfType );
[ Wt, ~, ~, ~ ] = phsWeights2d( Xt, Zt, xt, zt, [TxTop,TzTop], [NxTop,NzTop], ...
	rbforder, mn, np, K, rbfType );
[ halfWt, ~, ~, ~ ] = phsWeights2d( Xt, Zt, (xt+Xt(:,1))/2, (zt+Zt(:,1))/2, [TxTop,TzTop], [NxTop,NzTop], ...
	rbforder, mn, np, K, rbfType );
% [ Wb, Wt, halfWb, halfWt ] = getExtrapolationWeights( radb, radt, Xb, Zb, Xt, Zt, xb, zb, xt, zt, ...
	% n, sLayers, rbforder, mn, np, rbfType );
bigTx = bigTx(IDXb);
bigTz = bigTz(IDXb);
bigNx = bigNx(IDXb);
bigNz = bigNz(IDXb);

%interpolation and normal derivative on bottom boundary:
[ Wib, ~, Wndb, ~ ] = phsWeights2d( [xb,Xb], [zb,Zb], (xb+Xb(:,1))/2, (zb+Zb(:,1))/2, [Tx,Tz], [Nx,Nz], ...
	rbforder, mn, np, K, rbfType );

%interpolation and normal derivative on top boundary:
[ Wit, ~, Wndt, ~ ] = phsWeights2d( [xt,Xt], [zt,Zt], (xt+Xt(:,1))/2, (zt+Zt(:,1))/2, [TxTop,TzTop], [NxTop,NzTop], ...
	rbforder, mn, np, K, rbfType );

fprintf(1,'\ntime to generate weights = %g\n\n', toc );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Time-stepping:

if timeStep == 1
    
    tic
	
    if tPlot(1) ~= 0
		load( ['./matFiles/',saveName,'/',num2str(t(1)),'.mat'], 'U' )
    end
    
    mass = zeros( size(t) );
    initialMass = sum( 10^5*U(ind.m,1).^(Cv/Rd)./Rd./U(ind.m,4) .* cellAreas );
    mass(1) = initialMass;
    initialMassExact = sum( 10^5*U(ind.m,1).^(Cv/Rd)./Rd./U(ind.m,4) .* cellAreasExact );
    
    thetaMass = zeros( size(t) );
    initialThetaMass = sum( 10^5*U(ind.m,1).^(Cv/Rd)./Rd .* cellAreas );
    thetaMass(1) = initialThetaMass;
    initialThetaMassExact = sum( 10^5*U(ind.m,1).^(Cv/Rd)./Rd .* cellAreasExact );
    
    nul = zeros(size(U));

    for i = 0 : length(t)-1

        if i > 0
		
            if adamsBashforth == 1
                if i == 1
                    Fn2 = odeFun( t(i), U );
                    U = rk( t(i), U, dt, @odeFun, rkStages );
                elseif i == 2
                    Fn1 = odeFun( t(i), U );
                    U = rk( t(i), U, dt, @odeFun, rkStages );
                else
                    F0 = odeFun( t(i), U );
                    U = U + dt * ( 5/12*Fn2 - 4/3*Fn1 + 23/12*F0 );
                    Fn2 = Fn1;
                    Fn1 = F0;
                end
            elseif rbfMultistep == 1
                if i == 1
                    Fn4 = odeFun( t(i), U );
                    U = rk( t(i), U, dt, @odeFun, rkStages );
                elseif i == 2
                    Fn3 = odeFun( t(i), U );
                    U = rk( t(i), U, dt, @odeFun, rkStages );
                elseif i == 3
                    Fn2 = odeFun( t(i), U );
                    U = rk( t(i), U, dt, @odeFun, rkStages );
                elseif i == 4
                    Fn1 = odeFun( t(i), U );
                    U = rk( t(i), U, dt, @odeFun, rkStages );
                else
                    F0 = odeFun( t(i), U );
                    U = U + dt * ( 1/8*Fn4 - 5/24*Fn3 + 7/24*Fn2 - 23/24*Fn1 + 7/4*F0 );
                    Fn4 = Fn3;
                    Fn3 = Fn2;
                    Fn2 = Fn1;
                    Fn1 = F0;
                end
            else
                U = rk( t(i), U, dt, @odeFun, rkStages );
            end
			
            if useMassFixer==1
                %fix theta mass:
                rhoTheta = 10^5*U(ind.m,1).^(Cv/Rd)./Rd .* cellAreas;
                rhoTheta = rhoTheta + (initialThetaMass-sum(rhoTheta)).*areaRatios;
                rhoTheta = rhoTheta ./ cellAreas;
                U(ind.m,1) = (rhoTheta.*Rd/10^5).^(Rd/Cv);
                %fix regular mass:
                rho = 10^5*U(ind.m,1).^(Cv/Rd)./Rd./U(ind.m,4) .* cellAreas;
                rho = rho + (initialMass-sum(rho)).*areaRatios;
                rho = rho ./ cellAreas;
                U(ind.m,4) = 10^5*U(ind.m,1).^(Cv/Rd)./Rd./rho;
            end
			
            tmp = 10^5*U(ind.m,1).^(Cv/Rd)./Rd;
            mass(i+1) = sum( tmp./U(ind.m,4) .* cellAreas );
            thetaMass(i+1) = sum( tmp .* cellAreas );
			
        end

        if nnz(t(i+1)==tPlot)
            
            fprintf( 1, '\nt=%g, et=%g\nminPi=%g, maxPi=%g\nminU=%g, maxU=%g\nminW=%g, maxW=%g\nminTh=%g, maxTh=%g\ntotalMassChange=%g\ntotalThetaMassChange=%g\n', ...
                t(i+1), toc, min(U(ind.m,1)-piBar(ind.m)), max(U(ind.m,1)-piBar(ind.m)), ...
                min(U(ind.m,2)), max(U(ind.m,2)), ...
                min(U(ind.m,3)), max(U(ind.m,3)), ...
                min(U(ind.m,4)-thetaBar(ind.m)), max(U(ind.m,4)-thetaBar(ind.m)), ...
                ( mass(i+1) - initialMassExact ) / initialMassExact, ...
                ( thetaMass(i+1) - initialThetaMassExact ) / initialThetaMassExact );
            tic
            
            U = setGhostNodes( U );
            
            if saveResults == 1
                save( ['./matFiles/',saveName,'/',num2str(t(i+1)),'.mat'], 'U' )
            end
            
            plotContours( seeContours, testCase, xxc, zzc, xx, U, piBar, thetaBar, a, b, c, d, topoFunc, n, ind, Rd, Cp, Cv )

        end

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Y = odeFun( ~, U )

        %Assign values to ghost nodes:

        U = setGhostNodes( U );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Dx = Wx*U;
        Dz = Wz*U;
        HV = Whv*U;
%         HV = repmat(sqrt(U(ind.m,2).^2+U(ind.m,3).^2),1,4) .* (Whv*U);
%         V  = Wl*U(:,2:4);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Get RHS of system of ODEs for U = [pi,u,w,th]:

        U = U(ind.m,:);
        HV = repmat(-U(:,2),1,4).*Dx - repmat(U(:,3),1,4).*Dz + HV;
        HV(:,1:3) = HV(:,1:3) + [ -Rd/Cv*U(:,1).*(Dx(:,2)+Dz(:,3)), ...
            -Cp*U(:,4).*Dx(:,1), ...
            -Cp*U(:,4).*Dz(:,1)-g ];
        % HV(:,2:4) = HV(:,2:4) + V;
        
        Y = nul;
        Y(ind.m,:) = HV;

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function U = setGhostNodes( U )

        %periodic laterally in all variables (first time):
        U(ind.glr,:) = U(ind.rl,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%temporary variables:
		
		piTmp = U(:,1);
		uTmp  = U(:,2);
		wTmp  = U(:,3);
		thTmp = U(:,4);
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %theta on bottom (extrapolation, since no BC for theta):
        U(ind.gb,4) = sum( Wb.*thTmp(IDXb), 2 );
        %theta on top (extrapolation, since no BC for theta):
        U(ind.gt,4) = sum( Wt.*thTmp(IDXt), 2 );

        % %theta on bottom (imposing zero normal derivative in thetaPrime):
        % U(ind.gb,4) = ( U(ind.b,4) - thetaBar(ind.b) ) + thetaBar(ind.gb);
        % %theta on top:
        % U(ind.gt,4) = ( U(ind.t,4) - thetaBar(ind.t) ) + thetaBar(ind.gt);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%pi on bottom (enforcing momentum equation to high order):
		U(ind.gb,1) = ( -sum(Wndb(:,2:end).*piTmp(IDXb),2) - ...
			(Nz*g)./Cp./sum(halfWb.*thTmp(IDXb),2) ) ./ Wndb(:,1);
		%pi on top (enforcing momentum equation to high order):
		U(ind.gt,1) = ( -sum(Wndt(:,2:end).*piTmp(IDXt),2) - ...
			g./Cp./sum(halfWt.*thTmp(IDXt),2) ) ./ Wndt(:,1);

%         %pi on bottom (enforcing momentum equation):
%         U(ind.gb,1) = U(ind.b,1) + 2*alp.*(Nz*g)./Cp./sum(halfWb.*tmp(IDXb),2);
%         %pi on top (enforcing vertical momentum equation):
%         U(ind.gt,1) = U(ind.t,1) + 2*bet.*g./Cp./sum(halfWt.*tmp(IDXt),2);
		
% 		%pi on bottom (extrapolation):
% 		U(ind.gb,1) = sum( Wb.*piTmp(IDXb), 2 );
% 		%pi on top (extrapolation):
% 		U(ind.gt,1) = sum( Wt.*piTmp(IDXt), 2 );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %NOT using duT/dN=0 (correct free-slip condition):

        %uT on bottom:
        uT = uTmp(IDXb).*bigTx + wTmp(IDXb).*bigTz;
        uT = sum( Wb.*uT, 2 );                                %ghost values
		%uN on bottom:
        uN = uTmp(IDXb).*bigNx + wTmp(IDXb).*bigNz;
		uN = -sum(Wib(:,2:end).*uN,2) ./ Wib(:,1);			  %ghost values
		%combine to get u and w on bottom ghost nodes:
        U(ind.gb,2) = uT.*Tx + uN.*Nx;
        U(ind.gb,3) = uT.*Tz + uN.*Nz;
        %u on top:
        U(ind.gt,2) = sum( Wt.*uTmp(IDXt), 2 );
		%w on top:
		U(ind.gt,3) = -sum(Wit(:,2:end).*wTmp(IDXt),2) ./ Wit(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %periodic laterally in all variables (top/bottom ghost nodes only):
        U(ind.glr2,:) = U(ind.rl2,:);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
