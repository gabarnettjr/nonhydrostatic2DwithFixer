function [ IDX, rad, XC, ZC, h ] = ...
	getNearestNeighbors( simpleNeighbors, xxc, zzc, tmpX, tmpZ, nx, nz, n, sLayers, dx, dz )

if simpleNeighbors == 1
    
    x = xxc(:);
    z = zzc(:);
%     tmpX = xxc( (n+1)/2:end-(n-1)/2, 2:end-1 );
%     tmpZ = zzc( (n+1)/2:end-(n-1)/2, 2:end-1 );
    IDX = knnsearch( [x,(dx/dz)*z], [tmpX(:),(dx/dz)*tmpZ(:)], 'k', round(n*sLayers) );
	XC = x(IDX);
    ZC = z(IDX);
	rad = sqrt( (XC-repmat(tmpX(:),1,size(XC,2))).^2 + ...
		        (ZC-repmat(tmpZ(:),1,size(XC,2))).^2 );
    h = sum( rad(:,2:9), 2 ) / 8;
    rad = rad(:,end);
    
else

    %get nearest neighbors in 1D layers:
    [ idx, rad ] = knnsearch( xxc(:,2), xxc((n+1)/2:end-(n-1)/2,2), 'k', n );
    rad = repmat( rad(:,end), nz-1, 1 );
    rad = rad(:);

    if sLayers == 3
        %get three nearest layers:
        idx = [ idx, idx+(nx+n-2), idx+2*(nx+n-2) ];
        IDX = idx;
        for j = 1 : nz-2
            IDX = [ IDX; idx+j*(nx+n-2) ];
        end
    elseif sLayers == 5
        %get five nearest layers (unsymmetric):
        idx = [ idx, idx+(nx+n-2), idx+2*(nx+n-2), idx+3*(nx+n-2), idx+4*(nx+n-2) ];
        IDX = idx;
        IDX = [ IDX; idx ];
        for j = 1 : nz-4
            IDX = [ IDX; idx+j*(nx+n-2) ];
        end
        IDX = [ IDX; idx+(nz-4)*(nx+n-2) ];
    elseif sLayers == 7
        %get seven nearest layers (unsymmetric):
        idx = [ idx, idx+(nx+n-2), idx+2*(nx+n-2), idx+3*(nx+n-2), idx+4*(nx+n-2), idx+5*(nx+n-2), idx+6*(nx+n-2) ];
        IDX = idx;
        IDX = [ IDX; idx ];
        IDX = [ IDX; idx ];
        for j = 1 : nz-6
            IDX = [ IDX; idx+j*(nx+n-2) ];
        end
        IDX = [ IDX; idx+(nz-6)*(nx+n-2) ];
        IDX = [ IDX; idx+(nz-6)*(nx+n-2) ];
    end

    %write x and z centers in long vectors:
    xc = xxc(:);
    zc = zzc(:);

    %put stencils in rows:
    XC = xc(IDX);
    ZC = zc(IDX);
    
end