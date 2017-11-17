function ind = getIndexes( a, b, xxc, zzc, n )

x = xxc(:);
z = zzc(:);
tree = createns( [x,z] );

%main:
xx = xxc( (n+1)/2:end-(n-1)/2, 2:end-1 );
zz = zzc( (n+1)/2:end-(n-1)/2, 2:end-1 );
ind.m = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );

%left ghost nodes and corresponding right nodes:
xx = xxc( 1:(n-1)/2, 2:end-1 );
zz = zzc( 1:(n-1)/2, 2:end-1 );
ind.gl = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );
ind.r = knnsearch( tree, [x(ind.gl)+(b-a),z(ind.gl)], 'k', 1 );
%same, but top and bottom only:
xx = xxc( 1:(n-1)/2, [1,end] );
zz = zzc( 1:(n-1)/2, [1,end] );
ind.gl2 = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );
ind.r2 = knnsearch( tree, [x(ind.gl2)+(b-a),z(ind.gl2)], 'k', 1 );

%right ghost nodes and corresponding left nodes:
xx = xxc( end-(n-3)/2:end, 2:end-1 );
zz = zzc( end-(n-3)/2:end, 2:end-1 );
ind.gr = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );
ind.l = knnsearch( tree, [x(ind.gr)-(b-a),z(ind.gr)], 'k', 1 );
%same, but top and bottom only:
xx = xxc( end-(n-3)/2:end, [1,end] );
zz = zzc( end-(n-3)/2:end, [1,end] );
ind.gr2 = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );
ind.l2 = knnsearch( tree, [x(ind.gr2)-(b-a),z(ind.gr2)], 'k', 1 );

%bottom ghost nodes:
xx = xxc( (n+1)/2:end-(n-1)/2, 1 );
zz = zzc( (n+1)/2:end-(n-1)/2, 1 );
ind.gb = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );
%corresponding bottom nodes:
xx = xxc( (n+1)/2:end-(n-1)/2, 2 );
zz = zzc( (n+1)/2:end-(n-1)/2, 2 );
ind.b = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );

%top ghost nodes:
xx = xxc( (n+1)/2:end-(n-1)/2, end );
zz = zzc( (n+1)/2:end-(n-1)/2, end );
ind.gt = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );
%corresponding top nodes:
xx = xxc( (n+1)/2:end-(n-1)/2, end-1 );
zz = zzc( (n+1)/2:end-(n-1)/2, end-1 );
ind.t = knnsearch( tree, [xx(:),zz(:)], 'k', 1 );

ind.glr = [ ind.gl; ind.gr ];
ind.rl = [ ind.r; ind.l ];

ind.glr2 = [ ind.gl2; ind.gr2 ];
ind.rl2 = [ ind.r2; ind.l2 ];