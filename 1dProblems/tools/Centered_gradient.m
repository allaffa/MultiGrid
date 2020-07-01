function gradient = Centered_gradient(xh, uh)

    assert(size(xh,1)==size(uh,1));
    gradient = zeros(size(uh,1)-2,1);
    
    % We only compute the gradient for internal nodes. 
    % In case we have to refine, we will refine also on the adjacent area 
    % on the boundary

    for i = 2:size(xh,1)-1
        
        gradient(i-1) = ( uh(i+1) - uh(i-1) )/( xh(i+1) - xh(i-1) );
        
    end

    return;
    
end

