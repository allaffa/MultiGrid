function [omega, ratio] = compute_convergence_ratio( Hratiox, nx, Hratioy, ny )

    if( Hratiox > Hratioy )

        omega = 2 / ( sin( pi/(2 * Hratiox) ) * sin( pi/(2 * Hratiox) ) +...
        sin( pi/(2*nx) * (nx-1) ) * sin( pi/(2*nx) * (nx-1) ) + ...
        sin( pi/(2*ny) * (ny-1) ) * sin( pi/(2*ny) * (ny-1) ) );
    
    else
        
        omega = 2 / ( sin( pi/(2 * Hratioy) ) * sin( pi/(2 * Hratioy) ) +...
        sin( pi/(2*nx) * (nx-1) ) * sin( pi/(2*nx) * (nx-1) ) + ...
        sin( pi/(2*ny) * (ny-1) ) * sin( pi/(2*ny) * (ny-1) ) );
        
    end
    
    
    if (Hratiox > Hratioy)
    
        ratio = 1 - omega * sin( pi/(2 * Hratiox) ) * sin( pi/(2 * Hratiox) ) ;
        
    else
       
        ratio = 1 - omega * sin( pi/(2 * Hratioy) ) * sin( pi/(2 * Hratioy) ) ;
        
    end
    
return 
end