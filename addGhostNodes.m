function [ xxc, zzc, Nx, Nz, NxTop, NzTop, Tx, Tz, TxTop, TzTop, ...
	alp, bet, NC, bigTx, bigTz, bigNx, bigNz ] = addGhostNodes( xx, zz, xxc, zzc, a, b, n )

%get unit normal vectors on the bottom:
Nx = zz(1:end-1,1) - zz(2:end,1);
Nz = xx(2:end,1) - xx(1:end-1,1);
len = sqrt( Nx.^2 + Nz.^2 );
Nx = Nx ./ len;
Nz = Nz ./ len;

bigNx = [ Nx(end-(n-3)/2:end); Nx; Nx(1:(n-1)/2) ];
bigNx = repmat( bigNx, 1, size(xxc,2)+2 );
bigNx = bigNx(:);
bigNz = [ Nz(end-(n-3)/2:end); Nz; Nz(1:(n-1)/2) ];
bigNz = repmat( bigNz, 1, size(zzc,2)+2 );
bigNz = bigNz(:);

%get unit tangent vectors on the bottom:
Tx = xx(2:end,1) - xx(1:end-1,1);
Tz = zz(2:end,1) - zz(1:end-1,1);
len = sqrt( Tx.^2 + Tz.^2 );
Tx = Tx ./ len;
Tz = Tz ./ len;

bigTx = [ Tx(end-(n-3)/2:end); Tx; Tx(1:(n-1)/2) ];
bigTx = repmat( bigTx, 1, size(xxc,2)+2 );
bigTx = bigTx(:);
bigTz = [ Tz(end-(n-3)/2:end); Tz; Tz(1:(n-1)/2) ];
bigTz = repmat( bigTz, 1, size(zzc,2)+2 );
bigTz = bigTz(:);

%get bottom layer of ghost nodes:
x = xx(:,1);
z = zz(:,1);
xc = xxc(:,1) - x(1:end-1);
zc = zzc(:,1) - z(1:end-1);
alp = xc.*Nx + zc.*Nz;
% %ghost nodes lined up vertically:
% dz = zzc(:,2) - zzc(:,1);
% xxc = [ xxc(:,1), xxc ];
% zzc = [ zzc(:,1)-dz, zzc ];
%ghost nodes lined up in normal direction to boundary:
xxc = [ xxc(:,1)-2*alp.*Nx, xxc ];
zzc = [ zzc(:,1)-2*alp.*Nz, zzc ];

%get unit normal vectors on the top:
NxTop = zeros( size(xc) );
NzTop = ones( size(xc) );

%get unit tangent vectors on the top:
TxTop = ones( size(xc) );
TzTop = zeros( size(xc) );

%get top layer of ghost nodes:
x = xx(:,end);
z = zz(:,end);
xc = xxc(:,end) - x(1:end-1);
zc = zzc(:,end) - z(1:end-1);
bet = xc.*NxTop + zc.*NzTop;
% %ghost nodes lined up vertically:
% xxc = [ xxc, xxc(:,end) ];
% zzc = [ zzc, zzc(:,end)+dz ];
%ghost nodes lined up in normal direction to boundary:
xxc = [ xxc, xxc(:,end)-2*bet.*NxTop ];
zzc = [ zzc, zzc(:,end)-2*bet.*NzTop ];

%get left layers of ghost nodes:
xxc = [ xxc(end-(n-3)/2:end,:)-(b-a); xxc ];
zzc = [ zzc(end-(n-3)/2:end,:); zzc ];

%get right layers of ghost nodes:
xxc = [ xxc; xxc((n+1)/2:n-1,:)+(b-a) ];
zzc = [ zzc; zzc((n+1)/2:n-1,:) ];

%get total number of cells including ghost cells:
NC = length(xxc(:));