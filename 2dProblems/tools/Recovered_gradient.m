function Rg = Recovered_gradient(xh, yh, uh)

    assert( size(xh,1)*size(yh,1) == size(uh,1) );
    
    % We only compute the gradient for internal nodes. 
    % In case we have to refine, we will refine also on the adjacent area 
    % on the boundary
    
    nx = size(xh,1);
    ny = size(yh,1);
    
    Rg = zeros(nx*ny,1);

    for i = 2:nx-1
        
        for j = 2:ny-1
        
            kc = ( j - 1 ) * nx + i;
            ke = kc + 1;
            kw = kc - 1;
            kn = kc + nx;
            ks = kc - nx;
                  
            centered_fd_x = ( uh(ke) - uh(kw) )/( xh(i+1) - xh(i-1) );
            forward_fd_x = ( uh(ke) - uh(kc) )/( xh(i+1) - xh(i) );
            backward_fd_x = ( uh(kc) - uh(kw) )/( xh(i) - xh(i-1) );
            
            centered_fd_y = ( uh(kn) - uh(ks) )/( yh(j+1) - yh(j-1) );
            forward_fd_y = ( uh(kn) - uh(kc) )/( yh(j+1) - yh(j) );
            backward_fd_y = ( uh(kc) - uh(ks) )/( yh(j) - yh(j-1) );           

            Rg(kc, 1) = ( (centered_fd_x - forward_fd_x)^2 )*( xh(i+1) - xh(i) )+...
                       ( (centered_fd_x - backward_fd_x)^2 )*( xh(i) - xh(i-1) );
                   
            Rg(kc, 2) = ( (centered_fd_y - forward_fd_y)^2 )*( yh(j+1) - yh(j) )+...
                       ( (centered_fd_y - backward_fd_y)^2 )*( yh(j) - yh(j-1) );
                   

        end
        
    end

return 
end