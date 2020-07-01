function P = build_interpolant( Nx, nx, x1, x2, Ny, ny, y1, y2)

    hx = (x2 - x1)/(Nx - 1);   
    Hx = (x2 - x1)/(nx - 1);
    
    hy = (y2 - y1)/(Ny - 1);   
    Hy = (y2 - y1)/(ny - 1);
    
    N = Nx * Ny;
    n = nx * ny;
    
    ratiox = Hx/hx;
    ratioy = Hy/hy;
    
    P1 = sparse(N,n);
    
    % Interpolation for points on the edges of the coarse grid along
    % x-direction

    for i = 1:ny
        
        for j = 1 : nx-1

            for m = 0:ratiox-1

                P1( (i - 1) * (ratioy) * Nx + (j - 1) * ratiox + m + 1, (i - 1) * nx + j )   = ratiox - m;
                P1( (i - 1) * (ratioy) * Nx + (j - 1) * ratiox + m + 1, (i - 1) * nx + j+1 ) = m;

            end

        end
        
        P1( (i - 1) * (ratioy) * Nx + Nx, (i - 1) * nx + nx ) = ratiox;
        
    end
    
    P2 = sparse(N,n);
    
    % Interpolation for internal points along x-direction    
    
     for i = 1:ny-1
        
        for j = 1 : nx-1

            for m = 1:ratiox-1
                
                for k = 1:ratioy-1
                    
                    P2( (i - 1) * (ratioy) * Nx + k*Nx + (j - 1) * ratiox + m + 1, (i - 1) * nx + j )   = ratiox - m;
                    P2( (i - 1) * (ratioy) * Nx + k*Nx + (j - 1) * ratiox + m + 1, (i - 1) * nx + j+1 ) = m;
                    
                end

            end

        end
     end
        
     for i = 2:ny
        
        for j = 1 : nx-1

            for m = 1:ratiox-1
                
                for k = - (ratioy-1):-1
                    
                    P2( (i - 1) * (ratioy) * Nx + k*Nx + (j - 1) * ratiox + m + 1, (i - 1) * nx + j )   = ratiox - m;
                    P2( (i - 1) * (ratioy) * Nx + k*Nx + (j - 1) * ratiox + m + 1, (i - 1) * nx + j+1 ) = m;
                    
                end

            end

        end              
        
     end
    
     P3 = sparse(N,n);
    
    % Interpolation for points on  the edges of the coarse grid along
    % y-direction
     
    for j = 1:nx
        
        for i = 1:ny-1
            
            for m = 1:ratioy
                               
                P3( ( j - 1 ) * ratiox + 1 + ( i - 1 ) * ratioy * Nx + (m - 1) * Nx , (i - 1) * nx + j )        = ratioy - (m - 1) ;
                P3( ( j - 1 ) * ratiox + 1 + ( i - 1 ) * ratioy * Nx + (m - 1) * Nx , (i - 1) * nx + j + nx )   = m - 1;
                
            end
            
        end
        
        P3( ( j - 1 ) * ratiox + 1 + ( ny - 1 ) * ratioy * Nx , (ny - 1) * nx + j ) = ratioy;
        
    end
     
    % Interpolation for internal points along y-direction    
    
    P4 = sparse(N,n);
    
    for j = 1:nx-1

       for i = 1:ny-1

           for m = 2:ratioy
               
               for k = 1:ratiox-1

                   P4( ( j - 1 ) * ratiox + 1 + k + ( i - 1 ) * ratioy * Nx + (m - 1) * Nx , (i - 1) * nx + j )        = ratioy - (m - 1) ;
                   P4( ( j - 1 ) * ratiox + 1 + k + ( i - 1 ) * ratioy * Nx + (m - 1) * Nx , (i - 1) * nx + j + nx )   = m - 1;
               
               end

           end

       end  

    end
    
    for j = 2:nx

       for i = 1:ny-1

           for m = 2:ratioy
               
               for k = -(ratiox - 1):-1

                   P4( ( j - 1 ) * ratiox + 1 + k + ( i - 1 ) * ratioy * Nx + (m - 1) * Nx , (i - 1) * nx + j )        = ratioy - (m - 1) ;
                   P4( ( j - 1 ) * ratiox + 1 + k + ( i - 1 ) * ratioy * Nx + (m - 1) * Nx , (i - 1) * nx + j + nx )   = m - 1;
               
               end

           end

       end  

    end    
    
    Pb = ( ratioy * P1 + ratiox * P3 );
    
    Pi = P2 .* P4;
    
    P = 1/(ratiox * ratioy) * ( Pb + Pi );
    
    indices = P > 1.0;
    P(indices) = 1;
       
return
end