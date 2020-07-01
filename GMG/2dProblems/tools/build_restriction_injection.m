function R = build_restriction_injection( Nx, nx, x1, x2, Ny, ny, y1, y2 )

    hx = (x2 - x1)/(Nx - 1);   
    Hx = (x2 - x1)/(nx - 1);
    
    hy = (y2 - y1)/(Ny - 1);   
    Hy = (y2 - y1)/(ny - 1);
    
    N = Nx * Ny;
    n = nx * ny;
    
    R = sparse(n,N);
    
    ratiox = Hx/hx;
    ratioy = Hy/hy;
    
    countx = 0;
    county = 0;
    
    for i = 1:ny
        for j = 1:nx
            R(j + (i-1)*nx, countx*ratiox + county*Nx + 1) = 1;
            countx = countx + 1;
        end
        county = county + ratioy;
        countx = 0;
    end


return
end