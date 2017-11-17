function plotStencils( simpleNeighbors, seeStencils, xxc, zzc, Nx, Nz, NxTop, NzTop, alp, bet, XC, ZC, n, sLayers, nx, nz )

if seeStencils == 1

    plot( [xxc((n+1)/2:end-(n-1)/2,1),xxc((n+1)/2:end-(n-1)/2,1)+2*alp.*Nx].', ...
        [zzc((n+1)/2:end-(n-1)/2,1),zzc((n+1)/2:end-(n-1)/2,1)+2*alp.*Nz].', 'r-' )
    
    plot( [xxc((n+1)/2:end-(n-1)/2,end),xxc((n+1)/2:end-(n-1)/2,end)+2*bet.*NxTop].', ...
        [zzc((n+1)/2:end-(n-1)/2,end),zzc((n+1)/2:end-(n-1)/2,end)+2*bet.*NzTop].', 'r-' )
    
    tmp = xxc((n+1)/2:end-(n-1)/2,2:end-1);
    for i = 1 : length(tmp(:))
        plot( XC(i,:), ZC(i,:), 'ko' )
        if simpleNeighbors == 1
            plot( XC(i,1), ZC(i,1), 'ro' )
        else
            if sLayers == 3
                plot( XC(i,n+1), ZC(i,n+1), 'ro' )
            elseif sLayers == 5
                if i <= nx-1
                    plot( XC(i,n+1), ZC(i,n+1), 'ro' )
                elseif i > (nz-2)*(nx-1)
                    plot( XC(i,3*n+1), ZC(i,3*n+1), 'ro' )
                else
                    plot( XC(i,2*n+1), ZC(i,2*n+1), 'ro' )
                end
            elseif sLayers == 7
                if i <= nx-1
                    plot( XC(i,n+1), ZC(i,n+1), 'ro' )
                elseif i > nx-1 && i <=2*(nx-1)
                    plot( XC(i,2*n+1), ZC(i,2*n+1), 'ro' )
                elseif i > (nz-3)*(nx-1) && i <= (nz-2)*(nx-1)
                    plot( XC(i,4*n+1), ZC(i,4*n+1), 'ro' )
                elseif i > (nz-2)*(nx-1)
                    plot( XC(i,5*n+1), ZC(i,5*n+1), 'ro' )
                else
                    plot( XC(i,3*n+1), ZC(i,3*n+1), 'ro'  )
                end
            end
        end
        drawnow,pause
        plot( XC(i,:), ZC(i,:), 'wo' )
    end
    
end