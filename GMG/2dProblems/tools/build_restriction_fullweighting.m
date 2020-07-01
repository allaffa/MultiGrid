function R = build_restriction_fullweighting( Nx, nx, x1, x2, Ny, ny, y1, y2 )

    hx = (x2 - x1)/(Nx - 1);   
    Hx = (x2 - x1)/(nx - 1);
    
    hy = (y2 - y1)/(Ny - 1);   
    Hy = (y2 - y1)/(ny - 1);
    
    N = Nx * Ny;
    n = nx * ny;
    
    ratiox = Hx/hx;
    ratioy = Hy/hy;    
     
    denominator_x = ratiox * ratiox;
    denominator_y = ratioy * ratioy; 
    
    denominator = denominator_x * denominator_y;
    
    R1 = sparse(N,n);
    
    for j = 2:ny-1
    
        for i = 2:nx-1

            Fcenter = ( j-1 ) * ratioy * Nx + ( i-1 ) * ratiox + 1 ;
            Ccenter  = ( j-1 ) * nx + i;


            for k = - (ratioy - 1) : ratioy - 1
            
                for l = 1:(ratiox-1)
                    R1( Fcenter + k * Nx + l , Ccenter ) = ratiox-l;
                end

                R1( Fcenter + k * Nx , Ccenter ) = ratiox;
                
                for l = 1:(ratiox-1)
                    R1( Fcenter + k * Nx - l , Ccenter ) = ratiox-l;
                end
            
            end

        end
        
    end  
     
    
    R2 = sparse(N,n);
    
    for i = 2:nx-1
    
        for j = 2:ny-1

            Fcenter = ( j-1 ) * ratioy * Nx + ( i-1 ) * ratiox + 1 ;
            Ccenter  = i + ( j-1 ) * nx;

            for k = - (ratiox - 1) : ratiox - 1
            
                for l = 1:(ratioy-1)
                    R2( Fcenter + k + l * Nx , Ccenter ) = ratioy-l;
                end

                R2( Fcenter + k  , Ccenter ) = ratioy;
                
                for l = 1:(ratioy-1)
                    R2( Fcenter + k - l * Nx, Ccenter ) = ratioy-l;
                end
            
            end

        end
        
    end  
   
    
    Rtilde = R1 .* R2;
    
    Rtilde = 1/denominator * Rtilde;
    
    
    Rb = sparse(N,n);
    
    % Include the restriction coefficients for the boundary
    
    for i = 1 : nx
        
        % lower boundary
        Rb( (i-1)*ratiox + 1, i ) = 1;
        
        % upper boundary
        Rb( (ny-1)*ratioy*Nx + (i-1)*ratiox + 1 , (ny-1) * nx + i) = 1;
        
    end
    
    for j = 1 : ny
        
        % left boundary
        Rb( (j-1) * ratioy * Nx + 1 , (j-1)*nx + 1 ) = 1;
        
        % right boundary
        Rb( (j-1) * ratioy * Nx + Nx , j * nx ) = 1; 
        
    end

    R = (Rtilde + Rb)';
    
return
end