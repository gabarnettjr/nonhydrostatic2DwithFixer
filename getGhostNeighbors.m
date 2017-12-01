function [ IDXb, IDXt, Xb, Zb, Xt, Zt, xb, zb, xt, zt ] = ...
    getGhostNeighbors( xxc, zzc, n, sLayers, ind, dx, dz )
	
stencilSize = round( n * sLayers );

X = xxc(:);
Z = zzc(:);

x = xxc(:,2:end-1);
z = zzc(:,2:end-1);
x = x(:);
z = z(:);

xb = X(ind.gb);
zb = Z(ind.gb);
IDXb = knnsearch( [x,(dx/dz)*z], [xb,(dx/dz)*zb], 'k', stencilSize );
IDXb = IDXb.';
IDXb = knnsearch( [X,Z], [x(IDXb(:)),z(IDXb(:))], 'k', 1 );
IDXb = reshape( IDXb, stencilSize, length(xb) );
IDXb = IDXb.';
Xb = X(IDXb);
Zb = Z(IDXb);

xt = X(ind.gt);
zt = Z(ind.gt);
IDXt = knnsearch( [x,(dx/dz)*z], [xt,(dx/dz)*zt], 'k', stencilSize );
IDXt = IDXt.';
IDXt = knnsearch( [X,Z], [x(IDXt(:)),z(IDXt(:))], 'k', 1 );
IDXt = reshape( IDXt, stencilSize, length(xb) );
IDXt = IDXt.';
Xt = X(IDXt);
Zt = Z(IDXt);

% for i = 1 : size(Xb,1)
    % figure(10),clf
        % plot( X, Z, 'k.' )
        % hold( 'on' )
        % plot( Xb(i,:), Zb(i,:), 'ko' )
        % plot( Xb(i,1), Zb(i,1), 'ro' )
        % plot( xb(i), zb(i), 'bo' )
        % hold( 'off' )
        % axis( 'equal', 'tight', 'off' )
    % drawnow,pause
% end

% for i = 1 : size(Xt,1)
    % figure(10),clf
        % plot( X, Z, 'k.' )
        % hold( 'on' )
        % plot( Xt(i,:), Zt(i,:), 'ko' )
        % plot( Xt(i,1), Zt(i,1), 'ro' )
        % plot( xt(i), zt(i), 'bo' )
        % hold( 'off' )
        % axis( 'equal', 'tight', 'off' )
    % drawnow,pause
% end