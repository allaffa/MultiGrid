function [omega, ratio] = compute_convergence_ratio( Hratio, n )

    omega = 1 / ( sin( pi/(2 * Hratio) ) * sin( pi/(2 * Hratio) ) +...
        sin( pi/(2*n) * (n-1) ) * sin( pi/(2*n) * (n-1) ) );
    
    ratio = 2 * omega - 1;
    
return 
end
