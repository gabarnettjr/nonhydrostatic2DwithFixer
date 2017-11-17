function [ e1, e2 ] = getNewCoordinates( xxc, zzc, n, dx )

xxc = xxc( (n-1)/2:end-(n-3)/2, 2:end-1 );
zzc = zzc( (n-1)/2:end-(n-3)/2, 2:end-1 );

e11 = xxc( 3:end, : ) - xxc( 1:end-2, : );
e12 = zzc( 3:end, : ) - zzc( 1:end-2, : );

xxc = xxc(2:end-1,:);
zzc = zzc(2:end-1,:);

e1 = [ e11(:), e12(:) ];
e1 = e1 ./ repmat( sqrt(sum(e1.^2,2)), 1, 2 );
e2 = [ -e1(:,2), e1(:,1) ];

% figure(1),clf,hold('on'),plot(xxc(:),zzc(:),'k.','markerSize',10), ...
%     plot([xxc(:),xxc(:)+dx/2*e1(:,1)].',[zzc(:),zzc(:)+dx/2*e1(:,2)].','k'), ...
%     plot([xxc(:),xxc(:)+dx/2*e2(:,1)].',[zzc(:),zzc(:)+dx/2*e2(:,2)].','r'), ...
%     axis('equal')
% drawnow