function plotContours( seeContours, testCase, xxc, zzc, xx, U, piBar, thetaBar, a, b, c, d, topoFunc, n, ind, Rd, Cp, Cv )

if seeContours == 1
    
    xxc = xxc( (n+1)/2:end-(n-1)/2, 2:end-1 );
    zzc = zzc( (n+1)/2:end-(n-1)/2, 2:end-1 );
    U = U(ind.m,:);
    thetaBar = thetaBar(ind.m);
    piBar = piBar(ind.m);

    flag = 1;
    if strcmp( testCase, 'straka' ) || strcmp( testCase, 'strakaTopo' )
%         con = -30e-4 : 4e-4 : 30e-4;
%         con = -38 : 4 : 38;
%         con = -30 : 4 : 30;
        con = -20.5 : 1 : 10.5;
		cbar.location = 'north';
		cbar.color = 'black';
    elseif strcmp( testCase, 'movingStraka' ) || strcmp( testCase, 'movingStrakaTopo' )
        % con = -9e-3 : 2e-3 : 9e-3;
        % con = -27.5 : 5 : 67.5;
        % con = -30 : 4 : 30;
        con = -20.5 : 1 : 5.5;
		cbar.location = 'north';
		cbar.color = 'black';
    elseif strcmp( testCase, 'doubleStraka' ) || strcmp( testCase, 'doubleStrakaTopo' ) ...
            || strcmp( testCase, 'doubleStrakaSmooth' ) || strcmp( testCase, 'doubleStrakaTopoSmooth' )
        % con = -38e-4 : 4e-4 : 38e-4;
%         con = -42 : 4 : 42;
%         con = -34 : 4 : 34;
        con = -20.5 : 1 : 5.5;
		cbar.location = 'north';
		cbar.color = 'black';
    elseif strcmp( testCase, 'bubble' ) || strcmp( testCase, 'bubbleSmooth' ) ...
            || strcmp( testCase, 'bubbleTopo' ) || strcmp( testCase, 'bubbleTopoSmooth' )
        % con = -36e-5 : 3e-5 : 36e-5;
        % con = -9.5 : 1 : 9.5;
        % con = -11.5 : 1 : 11.5;
        con = -.1 : .2 : 1.9;
		cbar.location = 'east';
		cbar.color = 'white';
    elseif strcmp( testCase, 'igw' )
        con = (-2:.1:2)*1e-6;
        % con = (-.003:.0001:.003);
        % con = (-19:2:29)*1e-4;
        cbar.location = 'west';
        cbar.color = 'black';
    else
        flag = 0;
        nContours = 20;
    end
    cbar.fontSize = 30;

    figure(7),clf
        if flag == 1
            contourf( xxc, zzc, reshape(U(:,4)-thetaBar,size(xxc)), con, 'lineStyle', 'none' )
            colormap(parula(length(con)-1))
            caxis([min(con),max(con)])
            tmp = colorbar( cbar.location );
            set( tmp, 'color', cbar.color, 'fontSize', cbar.fontSize )
        else
            contourf( xxc, zzc, reshape(U(:,1)-piBar,size(xxc)), nContours, 'lineStyle', 'none' )
            colorbar
        end
        hold( 'on' )
        plot( xx(:,1), topoFunc(xx(:,1)), 'k' )
        plot( [b,b,a,a], [topoFunc(b),d,d,topoFunc(a)], 'k' )
        hold( 'off' )
        axis( 'tight', 'off' )
        axis( 'equal' )
    drawnow
    
end
