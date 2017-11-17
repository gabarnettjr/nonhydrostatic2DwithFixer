function plotComputationalDomain( seeDomain, xx, zz, xxc, zzc, nx, nz, dx, n, ind )

if seeDomain == 1

    lw = 1.5;
    ms = 8;

    figure(1),clf
    hold( 'on' )
    for j = 1 : nz+1
        plot( xxc(:,j), zzc(:,j), '.', 'lineWidth', lw, 'markerSize', ms )
    end
    for i = 1 : nx
        plot( xx(i,:), zz(i,:), 'k', 'lineWidth', lw )
    end
    for j = 1 : nz
        plot( xx(:,j), zz(:,j), 'k', 'lineWidth', lw )
    end
    axis( 'equal', 'tight', [xx(1,1)-(n+1)/2*dx,xx(end,end)+(n+1)/2*dx,zz(1,1)-(n+1)/2*dx,zz(end,end)+(n+1)/2*dx] )
    
%     x = xxc(:);
%     z = zzc(:);
%     plot( x,z,'k.', ...
%         x(ind.gr),z(ind.gr),'ro', ...
%         x(ind.l),z(ind.l),'rs', ...
%         x(ind.gl),z(ind.gl),'go', ...
%         x(ind.r),z(ind.r),'gs', ...
%         x(ind.gl2),z(ind.gl2),'bv', ...
%         x(ind.r2),z(ind.r2),'bh', ...
%         x(ind.m),z(ind.m),'k*'),axis('equal'),drawnow
    
end