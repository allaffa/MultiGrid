function Rg = Recovered_gradient(xh, uh)

    assert(size(xh,1)==size(uh,1));
    Rg = zeros(size(uh,1)-2,1);
    
    % We only compute the gradient for internal nodes. 
    % In case we have to refine, we will refine also on the adjacent area 
    % on the boundary

    for i = 2:size(xh,1)-1
        
        centered_fd = ( uh(i+1) - uh(i-1) )/( xh(i+1) - xh(i-1) );
        forward_fd = ( uh(i+1) - uh(i) )/( xh(i+1) - xh(i) );
        backward_fd = ( uh(i) - uh(i-1) )/( xh(i) - xh(i-1) );
        
        Rg(i-1) = ( (centered_fd - forward_fd)^2 )*( xh(i+1) - xh(i) )+...
                  ( (centered_fd - backward_fd)^2 )*( xh(i) - xh(i-1) );
        
    end

return 
end